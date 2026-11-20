"""Run the existing fixed provider worker with source-only imports and a pinned package path."""
import sys

if not (sys.flags.isolated and sys.flags.no_site and sys.dont_write_bytecode
        and sys.pycache_prefix == "/dev/null"):
    raise SystemExit("provider requires isolated source-only Python startup")
if len(sys.argv)<5 or sys.argv[1]!="--site" or not sys.argv[2].startswith("/"):
    raise SystemExit("provider requires the absolute prepared site-packages path")
site = sys.argv[2]
args = sys.argv[3:]
if args[:1] == ["-u"]:
    args=args[1:]
if len(args)<2 or args[0]!="-c":
    raise SystemExit("provider accepts the fixed worker command form")
sys.path.append(site)
sys.argv=["-c",*args[2:]]
exec(compile(args[1], "<prepared-provider-worker>", "exec"), {"__name__":"__main__"})
