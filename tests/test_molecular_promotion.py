"""Atomic promotion must stay within one bind mount, not just one st_dev."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]
definition=importlib.util.spec_from_file_location('molecular_native_promotion',ROOT/'scripts/pipeline/run_molecular_native.py')
worker=importlib.util.module_from_spec(definition);definition.loader.exec_module(worker)


def test_physical_alias_adapter_preserves_scientific_arguments(tmp_path):
    case=tmp_path/'case';(case/'dft').mkdir(parents=True);(case/'repairs').mkdir()
    install=worker.physical_installer(lambda *args:args,Path('/old'),case)
    args=install('/old/calcs/dft/id/neu','/old/results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs/id/neu/run/branch',
        '/old/results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs/id/neu/run','hash',
        '/old/results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs/id/neu/run/producer.py',0,1)
    assert args[0]==case/'dft/id/neu' and args[1]==case/'repairs/id/neu/run/branch'
    assert (args[3],args[5],args[6])==('hash',0,1)


@pytest.mark.parametrize('outside',['/archive/immutable','/old/calcs/dft/../outside'])
def test_undeclared_promotion_paths_are_rejected(tmp_path,outside):
    case=tmp_path/'case';case.mkdir()
    install=worker.physical_installer(lambda *args:args,Path('/old'),case)
    with pytest.raises(ValueError):install(outside,'/old/calcs/dft/id/neu','/old/calcs/dft/id','hash','/old/calcs/dft/id/producer',0,1)


def test_real_bind_mount_exdev_is_reproduced_and_fixed(tmp_path):
    case=tmp_path/'case';canonical=case/'dft/id/neu';branch=case/'repairs/run/branch'
    canonical.mkdir(parents=True);branch.mkdir(parents=True)
    (canonical/'marker').write_text('original');(branch/'marker').write_text('replacement')
    (case/'repairs/run/producer.py').write_text('source')
    prior=Path('/tmp/molecular_alias_control')
    code='''import errno,importlib.util,json,os,pathlib
spec=importlib.util.spec_from_file_location("tested_physical_promotion",WORKER);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
prior=pathlib.Path(PRIOR);case=pathlib.Path(CASE)
canonical=prior/"calcs/dft/id/neu";base=prior/"results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs/run";branch=base/"branch"
assert canonical.stat().st_dev==base.stat().st_dev
try:os.replace(canonical,base/"superseded")
except OSError as error:assert error.errno==errno.EXDEV
else:raise RuntimeError("The intended cross-mount failure was not reproduced")
def original(canonical,selected,base,digest,source,q,mult):
 assert (digest,q,mult)==("hash",0,1)
 os.replace(canonical,base/"superseded");os.replace(selected,canonical)
 return base/"superseded"
archive=m.physical_installer(original,prior,case)(canonical,branch,base,"hash",base/"producer.py",0,1)
assert (case/"dft/id/neu/marker").read_text()=="replacement"
assert (archive/"marker").read_text()=="original"
print(json.dumps(dict(original_error_errno=errno.EXDEV,same_st_dev=True,mapped_atomic_promotion_passed=True)))
'''
    program='WORKER='+repr(str(ROOT/'scripts/pipeline/run_molecular_native.py'))+'\nPRIOR='+repr(str(prior))+'\nCASE='+repr(str(case))+'\n'+code
    binary=ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap'
    argv=[str(binary),'--die-with-parent','--ro-bind','/','/','--tmpfs','/tmp','--bind',str(tmp_path),str(tmp_path),
        '--dir',str(prior/'calcs/dft'),'--dir',str(prior/'results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs'),
        '--bind',str(case/'dft'),str(prior/'calcs/dft'),'--bind',str(case/'repairs'),
        str(prior/'results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs'),
        '--proc','/proc','--',sys.executable,'-c',program]
    result=subprocess.run(argv,capture_output=True,text=True,timeout=30,check=True)
    assert json.loads(result.stdout)['mapped_atomic_promotion_passed']
