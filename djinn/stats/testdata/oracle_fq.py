"""Runs harpy's bx_stats_fq.py: usage: oracle_fq.py HARPY_UTILS PLATFORM INPUT.fq

Loads bx_stats_fq.py and harpy/common/parsers.py directly, so the rest of the harpy package is not needed."""
import importlib.util, os, sys, types

utils, platform, fastq = sys.argv[1], sys.argv[2], sys.argv[3]
for name in ("harpy", "harpy.common"):
    sys.modules[name] = types.ModuleType(name)
    sys.modules[name].__path__ = []

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

load("harpy.common.parsers", os.path.join(utils, "..", "common", "parsers.py"))
load("bx_stats_fq", os.path.join(utils, "bx_stats_fq.py")).bx_stats_fq(args=[platform, fastq], prog_name="bx_stats_fq")
