"""Boot-bound shared deadlines, with an append-only legacy-clock compatibility map."""
import ctypes
import datetime
from decimal import Decimal
from functools import lru_cache
from pathlib import Path
import time

from freeze import digest, dump, host_identity, load

NS = 1_000_000_000
DRIFT_NS = NS  # Operational clock guard, not a physical/statistical tolerance.
SCHEMA = 2
COMPANION = "clock-compatibility-v1.json"


class Budget(Exception):
    pass


@lru_cache(maxsize=1)
def _mach():
    lib = ctypes.CDLL("/usr/lib/libSystem.B.dylib", use_errno=True)
    lib.mach_continuous_time.restype = ctypes.c_uint64
    lib.mach_continuous_time.argtypes = []
    lib.sysctlbyname.argtypes = [ctypes.c_char_p, ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p, ctypes.c_size_t]
    lib.sysctlbyname.restype = ctypes.c_int
    return lib


def _sysctl(name, size):
    value = ctypes.create_string_buffer(size)
    count = ctypes.c_size_t(size)
    if _mach().sysctlbyname(name.encode(), value, ctypes.byref(count), None, 0):
        raise Budget("host clock identity unavailable: " + name + " errno=" + str(ctypes.get_errno()))
    return value.raw[:count.value]


@lru_cache(maxsize=1)
def clock_identity():
    # A live process cannot survive a reboot. Cache only its boot/source identity,
    # never its advancing clock samples. A new process reads the boot UUID again.
    class Timebase(ctypes.Structure):
        _fields_ = [("numer", ctypes.c_uint32), ("denom", ctypes.c_uint32)]
    tb = Timebase()
    if _mach().mach_timebase_info(ctypes.byref(tb)) or not tb.numer or not tb.denom:
        raise Budget("Mach timebase unavailable")
    boot = _sysctl("kern.bootsessionuuid", 128).rstrip(b"\0").decode("ascii")
    if not boot:
        raise Budget("invalid host boot identity")
    return dict(schema=SCHEMA, api="mach_continuous_time", library="/usr/lib/libSystem.B.dylib",
                numer=tb.numer, denom=tb.denom, boot_uuid=boot,
                host=host_identity(),
                implementation_sha256=digest(__file__), drift_guard_ns=DRIFT_NS)


def boot_utc_ns():
    class Timeval(ctypes.Structure):
        _fields_ = [("sec", ctypes.c_long), ("usec", ctypes.c_int)]
    tv = Timeval.from_buffer_copy(_sysctl("kern.boottime", ctypes.sizeof(Timeval)))
    if tv.sec <= 0 or not 0 <= tv.usec < 1_000_000:
        raise Budget("invalid boot UTC observation")
    return tv.sec*NS+tv.usec*1000


def shared_ns():
    identity = clock_identity()
    return _mach().mach_continuous_time()*identity["numer"]//identity["denom"]


def sample():
    identity = clock_identity()
    boot = boot_utc_ns()
    before = shared_ns()
    wall = time.time_ns()
    after = shared_ns()
    if after < before:
        raise Budget("shared clock reversed while sampling")
    return dict(before_ns=before, utc_ns=wall, after_ns=after, identity=identity, boot_utc_ns=boot)


def _utc_ns(value):
    # Floor the exact retained binary float; conversion cannot extend its deadline.
    return int(Decimal.from_float(float(value))*NS)


def new_clock(seconds=7200, observed=None):
    ref = sample() if observed is None else observed
    duration = int(seconds*NS)
    return dict(schema=SCHEMA, started_utc=datetime.datetime.fromtimestamp(
                    ref["utc_ns"]/NS, datetime.timezone.utc).isoformat(),
                origin_start_utc_ns=ref["utc_ns"], reference=ref,
                deadline_utc=(ref["utc_ns"]+duration)/NS,
                deadline_utc_ns=ref["utc_ns"]+duration,
                deadline_shared_ns=ref["before_ns"]+duration,
                elapsed_basis="boot-shared", legacy_monotonic_elapsed_seconds=None)


def _check_sample(clock, now):
    if clock.get("schema") != SCHEMA:
        raise Budget("legacy process-relative clock requires reviewed compatibility mapping")
    ref = clock["reference"]
    if now["identity"] != ref["identity"]:
        raise Budget("host clock boot/source identity changed")
    if (now["before_ns"] < ref["before_ns"] or now["after_ns"] < now["before_ns"]
            or now["utc_ns"] < ref["utc_ns"]):
        raise Budget("host clock discontinuity")
    elapsed_wall = now["utc_ns"]-ref["utc_ns"]
    lower = now["before_ns"]-ref["after_ns"]-DRIFT_NS
    upper = now["after_ns"]-ref["before_ns"]+DRIFT_NS
    if not lower <= elapsed_wall <= upper:
        raise Budget("UTC/shared clock divergence exceeds frozen guard")


def remaining(clock, observed=None):
    now = sample() if observed is None else observed
    _check_sample(clock, now)
    return min(clock["deadline_utc_ns"]-now["utc_ns"],
               clock["deadline_shared_ns"]-now["after_ns"])/NS


def limited_clock(clock, seconds):
    now = sample()
    _check_sample(clock, now)
    duration = max(0, int(seconds*NS))
    return dict(clock, deadline_utc_ns=min(clock["deadline_utc_ns"],now["utc_ns"]+duration),
                deadline_shared_ns=min(clock["deadline_shared_ns"],now["before_ns"]+duration),
                deadline_utc=min(clock["deadline_utc"],(now["utc_ns"]+duration)/NS))


def elapsed(clock, observed=None):
    now = sample() if observed is None else observed
    try:
        _check_sample(clock, now)
    except Budget as exc:
        return dict(seconds=None,basis="unknown",reason=str(exc),legacy_monotonic_elapsed_seconds=None)
    shared = (now["before_ns"]-clock["reference"]["before_ns"])/NS
    utc_elapsed = (now["utc_ns"]-clock["origin_start_utc_ns"])/NS
    migrated = clock["elapsed_basis"] == "legacy-prefix-UTC-derived"
    return dict(seconds=utc_elapsed if migrated else shared,basis=clock["elapsed_basis"],
                shared_seconds_since_reference=shared,utc_derived_seconds_since_origin=utc_elapsed,
                legacy_monotonic_elapsed_seconds=None,legacy_boot_identity=None,
                verified_boot_uuid=clock["reference"]["identity"].get("boot_uuid"))


def legacy_snapshot(ledger):
    """Read-only legacy binding; used in preparation and checked before migration."""
    anchor = ledger/"absolute-clock.json"
    if not anchor.exists():
        return None
    old = load(anchor)
    if old.get("schema") == SCHEMA:
        return None
    required = {"started_utc", "start_monotonic", "deadline_monotonic", "deadline_utc"}
    if set(old) != required:
        raise Budget("unknown legacy clock schema")
    files = {str(anchor):digest(anchor)}
    for path in sorted(ledger.glob("epoch-*.json")):
        entry = load(path)
        prior = Path(entry["out"])
        for item in [path, prior/"clock.json", prior/"protocol.json", prior/"external-epoch.json"]:
            files[str(item)] = digest(item)
        clock = load(prior/"clock.json")
        if any(clock.get(k) != v for k,v in old.items()):
            raise Budget("legacy epoch clock differs from original anchor")
        if digest(prior/"protocol.json") != entry["internal"] or digest(prior/"external-epoch.json") != entry["external"]:
            raise Budget("legacy epoch registration changed")
        for folder in (prior/"jobs").iterdir():
            if not folder.is_dir() or set(p.name for p in folder.iterdir()) != {"denied.json"}:
                raise Budget("legacy clock migration requires no numerical launch")
            denial = load(folder/"denied.json")
            if denial.get("launched") is not False:
                raise Budget("legacy numerical work cannot be reinterpreted")
            files[str(folder/"denied.json")] = digest(folder/"denied.json")
    return dict(anchor_sha256=digest(anchor),anchor=old,files=files)


def mapped_clock(legacy, observed):
    old = legacy["anchor"]
    deadline = _utc_ns(old["deadline_utc"])
    origin = int(Decimal(str(datetime.datetime.fromisoformat(old["started_utc"]).timestamp()))*NS)
    if observed["utc_ns"] >= deadline:
        raise Budget("original absolute deadline expired; cannot restart")
    if observed["utc_ns"] < origin or observed["boot_utc_ns"] > origin:
        raise Budget("legacy UTC discontinuity or reboot after original start")
    clock = new_clock(observed=observed)
    clock.update(started_utc=old["started_utc"],origin_start_utc_ns=origin,
                 deadline_utc=old["deadline_utc"],deadline_utc_ns=deadline,
                 deadline_shared_ns=observed["before_ns"]+deadline-observed["utc_ns"],
                 elapsed_basis="legacy-prefix-UTC-derived",legacy_boot_identity=None)
    # v5 did not record a boot identity. An earlier kernel boot timestamp is a
    # rejection aid only, not retrospective proof of the unknown legacy prefix.
    return clock


def verify_mapping(mapping, legacy):
    if mapping.get("legacy") != legacy or mapping.get("clock") != mapped_clock(legacy,mapping["observed"]):
        raise Budget("clock compatibility mapping changed")
    if mapping["observed"]["identity"] != clock_identity():
        raise Budget("mapped clock boot/source identity changed")
    if remaining(mapping["clock"],mapping["adopted_observed"]) <= 0:
        raise Budget("original absolute deadline expired at mapping adoption")
    return mapping["clock"]


def verify_legacy_files(legacy):
    if legacy is not None:
        for name,sha in legacy["files"].items():
            if digest(name) != sha:
                raise Budget("legacy clock/registration/denial changed")


def active_clock(ledger):
    anchor = ledger/"absolute-clock.json"
    original = load(anchor)
    if original.get("schema") == SCHEMA:
        return dict(original,absolute_anchor_sha256=digest(anchor))
    companion = ledger/COMPANION
    if not companion.exists():
        raise Budget("legacy process-relative clock requires reviewed compatibility mapping")
    mapping = load(companion)
    # Existing registrations are append-only after migration; recheck all original bytes.
    legacy = mapping["legacy"]
    if legacy["anchor_sha256"] != digest(anchor) or legacy["anchor"] != original:
        raise Budget("original clock anchor changed")
    verify_legacy_files(legacy)
    clock = verify_mapping(mapping,legacy)
    return dict(clock,absolute_anchor_sha256=digest(anchor),compatibility_sha256=digest(companion))


def migrate(ledger, expected_legacy, expected_identity, prepared_observation):
    """Explicit, one-time operator action AFTER source review. Never called by preparation."""
    if (ledger/COMPANION).exists():
        raise ValueError("clock compatibility mapping already exists")
    if expected_legacy is None or legacy_snapshot(ledger) != expected_legacy:
        raise Budget("legacy anchor/registrations differ from prepared compatibility binding")
    if prepared_observation["identity"] != expected_identity:
        raise Budget("prepared clock boot/source identity changed")
    clock = mapped_clock(expected_legacy,prepared_observation)
    adopted = sample()
    if remaining(clock,adopted) <= 0:
        raise Budget("original absolute deadline expired; cannot restart")
    mapping = dict(legacy=expected_legacy,observed=prepared_observation,
                   adopted_observed=adopted,clock=clock,
                   policy="original UTC deadline and conservative shared-before deadline; no reset")
    dump(ledger/COMPANION,mapping)
    return mapping


def legacy_elapsed(original, wall_ns=None):
    """Reporting only. This UTC-derived prefix never authorizes admission or migration."""
    now = time.time_ns() if wall_ns is None else wall_ns
    origin = int(Decimal(str(datetime.datetime.fromisoformat(original["started_utc"]).timestamp()))*NS)
    return dict(seconds=max(0,now-origin)/NS,basis="legacy-prefix-UTC-derived; unmapped",
                legacy_monotonic_elapsed_seconds=None,legacy_boot_identity=None,
                original_deadline_utc=original["deadline_utc"],
                original_deadline_expired=now>=_utc_ns(original["deadline_utc"]),
                admission="blocked: compatibility mapping absent")
