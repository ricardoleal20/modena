"""Independent-process clock controls. All anchors/epochs are temporary and synthetic."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
LAUNCH=ROOT/'bench/b07_validation/launch.sh'


class IndependentClock(unittest.TestCase):
    def run_python(self, code, *args):
        result=subprocess.run(['/bin/sh',str(LAUNCH),sys.executable,'-c',
            'import sys;sys.path.insert(0,'+repr(str(ROOT/'bench/b07_validation'))+');'+code,*map(str,args)],
            capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        return json.loads(result.stdout)

    def test_fresh_clock_survives_independent_process_admission(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td).resolve(); out=root/'epoch'; out.mkdir(); (out/'jobs').mkdir(); (out/'inputs').mkdir()
            (out/'protocol.json').write_text('{}'); (out/'external-epoch.json').write_text('{}')
            code="""import campaign as c,json,time
from pathlib import Path
from unittest.mock import patch
out=Path(sys.argv[1]);c.LEDGER=out.parent/'ledger'
if hasattr(c,'clock_identity'):
 (out/'inputs/settings.json').write_text(json.dumps({'clock_identity':c.clock_identity()}))
time.sleep(.25)
with patch.object(c,'verify',return_value=({},{})),patch.object(c,'jobs',return_value={}): c.start(out)
print(json.dumps(c.load(out/'clock.json')))
"""
            first=self.run_python(code,out)
            code="""import campaign as c,json
from pathlib import Path
try: value=dict(status='admitted',remaining=c.remaining(c.load(Path(sys.argv[1])/'clock.json')))
except c.Budget as exc: value=dict(status='denied',reason=str(exc))
print(json.dumps(value))
"""
            second=self.run_python(code,out)
            self.assertEqual(second['status'],'admitted',second)
            self.assertGreater(second['remaining'],0)
            self.assertLessEqual(second['remaining'],7200)


import datetime
import copy
from unittest.mock import patch
import campaign
import timing
from freeze import digest,dump,load


class ClockMapping(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name).resolve(); self.ledger=self.root/'ledger'; self.ledger.mkdir()
        self.prior=self.root/'old'; self.prior.mkdir(); (self.prior/'jobs/deterministic').mkdir(parents=True)
        observed=timing.sample(); start=observed['utc_ns']/timing.NS-30
        self.old=dict(started_utc=datetime.datetime.fromtimestamp(start,datetime.timezone.utc).isoformat(),
            start_monotonic=2.8,deadline_monotonic=122.8,deadline_utc=start+120)
        dump(self.ledger/'absolute-clock.json',self.old)
        dump(self.prior/'protocol.json',{'synthetic':'old internal'})
        dump(self.prior/'external-epoch.json',{'synthetic':'old external'})
        pins=dict(internal=digest(self.prior/'protocol.json'),external=digest(self.prior/'external-epoch.json'))
        dump(self.prior/'clock.json',dict(self.old,**pins))
        dump(self.ledger/'epoch-original.json',dict(out=str(self.prior),**pins))
        dump(self.prior/'jobs/deterministic/denied.json',dict(name='deterministic',status='BudgetExceeded',launched=False,reason='synthetic old clock discontinuity'))
        self.expected=timing.legacy_snapshot(self.ledger)
        self.prepared=timing.sample()
        self.before={p:Path(p).read_bytes() for p in self.expected['files']}

    def assert_preserved(self):
        self.assertEqual(self.before,{p:Path(p).read_bytes() for p in self.before})

    def test_mapping_is_conservative_one_time_and_preserves_original_bytes(self):
        mapping=timing.migrate(self.ledger,self.expected,timing.clock_identity(),self.prepared)
        obs=mapping['observed']; clock=mapping['clock']
        self.assertEqual(clock['deadline_utc'],self.old['deadline_utc'])
        self.assertEqual(clock['deadline_shared_ns'],obs['before_ns']+clock['deadline_utc_ns']-obs['utc_ns'])
        self.assertLessEqual(clock['deadline_shared_ns'],obs['after_ns']+clock['deadline_utc_ns']-obs['utc_ns'])
        self.assertGreater(timing.remaining(timing.active_clock(self.ledger)),0)
        details=timing.elapsed(clock)
        self.assertEqual(details['basis'],'legacy-prefix-UTC-derived')
        self.assertIsNone(details['legacy_monotonic_elapsed_seconds'])
        self.assertIsNone(clock['legacy_boot_identity'])
        self.assertEqual(mapping['observed'],self.prepared)
        self.assertGreater(details['utc_derived_seconds_since_origin'],30)
        with self.assertRaises(ValueError): timing.migrate(self.ledger,self.expected,timing.clock_identity(),self.prepared)
        self.assert_preserved()

    def test_review_interval_cannot_recalibrate_a_clock_jump(self):
        now=dict(self.prepared,
            before_ns=self.prepared['before_ns']+10*timing.NS,
            after_ns=self.prepared['after_ns']+10*timing.NS,
            utc_ns=self.prepared['utc_ns']+7*timing.NS)
        with patch.object(timing,'sample',return_value=now),self.assertRaisesRegex(timing.Budget,'divergence'):
            timing.migrate(self.ledger,self.expected,timing.clock_identity(),self.prepared)
        self.assertFalse((self.ledger/timing.COMPANION).exists()); self.assert_preserved()

    def test_admission_rejects_self_consistent_recalibrated_mapping(self):
        mapping=timing.migrate(self.ledger,self.expected,timing.clock_identity(),self.prepared)
        out=self.root/'new'; (out/'inputs').mkdir(parents=True)
        dump(out/'inputs/settings.json',dict(legacy_clock=self.expected,clock_observation=self.prepared))
        with patch.object(campaign,'LEDGER',self.ledger):
            campaign.prepared_clock(out)
            mapping['observed']=dict(self.prepared,utc_ns=self.prepared['utc_ns']-1)
            mapping['clock']=timing.mapped_clock(self.expected,mapping['observed'])
            (self.ledger/timing.COMPANION).write_text(json.dumps(mapping))
            with self.assertRaisesRegex(timing.Budget,'prepared legacy binding'):
                campaign.prepared_clock(out)
        self.assert_preserved()

    def test_mapped_clock_survives_independent_processes(self):
        helper=IndependentClock('runTest')
        code="import timing,json;from pathlib import Path;p=Path(sys.argv[1]);m=timing.migrate(p,timing.legacy_snapshot(p),timing.clock_identity(),timing.sample());print(json.dumps(m['clock']))"
        helper.run_python(code,self.ledger)
        code="import timing,json;from pathlib import Path;c=timing.active_clock(Path(sys.argv[1]));print(json.dumps({'remaining':timing.remaining(c),'elapsed':timing.elapsed(c)}))"
        actual=helper.run_python(code,self.ledger)
        self.assertGreater(actual['remaining'],0)
        self.assertLess(actual['remaining'],90)
        self.assertIsNone(actual['elapsed']['legacy_monotonic_elapsed_seconds'])
        self.assert_preserved()

    def test_expired_legacy_cannot_map_or_restart(self):
        now=timing.sample(); delta=timing._utc_ns(self.old['deadline_utc'])+1-now['utc_ns']
        now={k:v+delta if k in ('before_ns','after_ns','utc_ns') else v for k,v in now.items()}
        with patch.object(timing,'sample',return_value=now),self.assertRaisesRegex(timing.Budget,'expired'):
            timing.migrate(self.ledger,self.expected,timing.clock_identity(),self.prepared)
        self.assertFalse((self.ledger/timing.COMPANION).exists())
        out=self.root/'new'; out.mkdir(); (out/'jobs').mkdir(); (out/'inputs').mkdir()
        dump(out/'inputs/settings.json',{'clock_identity':timing.clock_identity()})
        with patch.object(campaign,'LEDGER',self.ledger),patch.object(campaign,'verify',return_value=({},{})),patch.object(campaign,'jobs',return_value={}):
            with self.assertRaises(timing.Budget): campaign.start(out)
        self.assertFalse((out/'clock.json').exists())
        self.assertEqual(len(list(self.ledger.glob('epoch-*.json'))),1)
        self.assert_preserved()

    def test_tampered_mapping_or_original_registration_fails_closed(self):
        timing.migrate(self.ledger,self.expected,timing.clock_identity(),self.prepared)
        file=self.ledger/timing.COMPANION; original=file.read_bytes(); bad=load(file)
        bad['clock']['deadline_shared_ns']+=1
        file.write_text(json.dumps(bad))
        with self.assertRaisesRegex(timing.Budget,'mapping changed'): timing.active_clock(self.ledger)
        file.write_bytes(original)
        reg=self.ledger/'epoch-original.json'; saved=reg.read_bytes(); reg.write_text(saved.decode()+' ')
        with self.assertRaisesRegex(timing.Budget,'registration'): timing.active_clock(self.ledger)
        reg.write_bytes(saved); self.assert_preserved()

    def test_changed_prepared_anchor_and_legacy_reboot_reject_without_mapping(self):
        bad=copy.deepcopy(self.expected); bad['anchor_sha256']='0'*64
        with self.assertRaisesRegex(timing.Budget,'prepared'): timing.migrate(self.ledger,bad,timing.clock_identity(),self.prepared)
        now=timing.sample(); now['boot_utc_ns']=now['utc_ns']
        with patch.object(timing,'sample',return_value=now),self.assertRaisesRegex(timing.Budget,'reboot'):
            timing.migrate(self.ledger,self.expected,timing.clock_identity(),now)
        self.assertFalse((self.ledger/timing.COMPANION).exists()); self.assert_preserved()

    def test_final_report_retains_original_denial_without_mapping(self):
        out=self.root/'new'; out.mkdir(); (out/'inputs').mkdir(); (out/'jobs').mkdir()
        dump(out/'inputs/settings.json',{'prior_projected_draws':[]})
        external={'expected_gate_families':{'synthetic':['pending-physical']}}
        with patch.object(campaign,'LEDGER',self.ledger),patch.object(campaign,'verify',return_value=({},external)),patch.object(campaign,'jobs',return_value={}):
            campaign.report(out)
        report=max((out/'reports').iterdir()); ledger=load(report/'gate-ledger.json'); work=load(report/'work.json')
        self.assertEqual(ledger['prior_epoch_denials'][0]['denial']['status'],'BudgetExceeded')
        self.assertFalse(work['prior_epoch_denials'][0]['denial']['launched'])
        self.assertIsNone(work['combined_elapsed_details']['legacy_monotonic_elapsed_seconds'])
        self.assertIn('UTC-derived',work['combined_elapsed_details']['basis'])
        self.assertEqual(work['cumulative_work']['all_molecule_steps'],0)
        self.assertFalse((self.ledger/timing.COMPANION).exists()); self.assert_preserved()


class ClockBoundaries(unittest.TestCase):
    def setUp(self):
        self.ref=dict(before_ns=100*timing.NS,after_ns=100*timing.NS+20,
                      utc_ns=1_700_000_000*timing.NS,identity={'boot':'synthetic'},boot_utc_ns=1)
        self.clock=timing.new_clock(100,self.ref)

    def test_frozen_guard_boundaries_and_utc_jumps(self):
        now=dict(self.ref,before_ns=self.ref['before_ns']+5*timing.NS,after_ns=self.ref['after_ns']+5*timing.NS)
        for delta in [now['before_ns']-self.ref['after_ns']-timing.DRIFT_NS,
                      now['after_ns']-self.ref['before_ns']+timing.DRIFT_NS]:
            self.assertGreater(timing.remaining(self.clock,dict(now,utc_ns=self.ref['utc_ns']+delta)),0)
        for delta in [5*timing.NS-timing.DRIFT_NS-21,5*timing.NS+timing.DRIFT_NS+21]:
            with self.assertRaisesRegex(timing.Budget,'divergence'):
                timing.remaining(self.clock,dict(now,utc_ns=self.ref['utc_ns']+delta))
        with self.assertRaisesRegex(timing.Budget,'discontinuity'):
            timing.remaining(self.clock,dict(now,utc_ns=self.ref['utc_ns']-1))

    def test_shared_reversal_boot_change_and_expiry(self):
        with self.assertRaisesRegex(timing.Budget,'discontinuity'):
            timing.remaining(self.clock,dict(self.ref,before_ns=self.ref['before_ns']-1))
        with self.assertRaisesRegex(timing.Budget,'identity'):
            timing.remaining(self.clock,dict(self.ref,identity={'boot':'rebooted'}))
        now={k:v+101*timing.NS if k in ('before_ns','after_ns','utc_ns') else v for k,v in self.ref.items()}
        self.assertLess(timing.remaining(self.clock,now),0)

    def test_clock_query_failure_is_visible(self):
        timing.clock_identity.cache_clear()
        try:
            with patch.object(timing,'_sysctl',side_effect=timing.Budget('unavailable')):
                with self.assertRaisesRegex(timing.Budget,'unavailable'): timing.clock_identity()
        finally:
            timing.clock_identity.cache_clear()

    def test_local_duration_is_a_difference_within_one_process(self):
        result=IndependentClock('runTest').run_python("import time,json;a=time.monotonic();time.sleep(.1);print(json.dumps({'elapsed':time.monotonic()-a}))")
        self.assertGreaterEqual(result['elapsed'],.09)
        self.assertLess(result['elapsed'],1)
