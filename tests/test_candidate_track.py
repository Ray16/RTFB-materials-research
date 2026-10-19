"""Candidate preparation contracts: no duplicate registration or inherited parent metrics."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest
from rdkit import Chem
from redox.build.contracts import validate_group

ROOT=Path(__file__).resolve().parents[1]


def _prepare_module(monkeypatch,groups,existing=()):
    spec=importlib.util.spec_from_file_location('candidate_preparation_test',ROOT/'scripts/mining/prepare_next_batch.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    monkeypatch.setattr(mod,'SELECTED',[g['id'] for g in groups])
    evidence={g['id']:dict(type='proposal',route_source='https://example.org/source',
       route_status='proposed',route_note='not experimentally validated',direct_merrifield_route_verified=False,
       experimental_parent_smiles='Cn1c(=O)c2ccccc2c1=O',experimental_parent_reference='parent DOI',
       parent_E_V=-1.8,parent_lambda_eV=.1,parent_solubility_M=1.5) for g in groups}
    monkeypatch.setattr(mod,'EVIDENCE',evidence)
    monkeypatch.setattr(mod,'read_manifest',lambda:[])
    monkeypatch.setattr(mod,'load_config',lambda name:SimpleNamespace(FARADAY=96485.33212,GROUPS=list(existing)))
    return mod


def _group(gid='new_imide'):
    return dict(id=gid,name='N-benzyl phthalimide',family='mono-imide',
                frag='[*:1]N1C(=O)c2ccccc2C1=O',states=[['neu',0,1,0],['red',-1,2,-1]])


def test_preparation_excludes_registered_graph_under_an_alias(tmp_path,monkeypatch):
    registered=_group('already_registered');new=_group('new_alias')
    mod=_prepare_module(monkeypatch,[new],existing=[registered])
    inp=tmp_path/'curated.json';inp.write_text(json.dumps([new]));out=tmp_path/'output'
    monkeypatch.setattr(sys,'argv',['prepare','--curated',str(inp),'--out',str(out)])
    with pytest.raises(SystemExit) as exc:mod.main()
    assert exc.value.code==2
    assert not out.exists()


def test_preparation_preserves_witness_but_never_inherits_properties(tmp_path,monkeypatch):
    g=_group();mod=_prepare_module(monkeypatch,[g])
    inp=tmp_path/'curated.json';inp.write_text(json.dumps([g]));out=tmp_path/'output'
    monkeypatch.setattr(sys,'argv',['prepare','--curated',str(inp),'--out',str(out)])
    mod.main()
    import csv
    row=next(csv.DictReader(open(out/'screening_inputs.csv')))
    assert row['experimental_parent_smiles']
    assert 'E_V' not in row and 'lambda_eV' not in row and 'solubility_M' not in row
    assert all('parent_E' not in key and 'parent_lambda' not in key for key in row)
    config=json.loads((out/'batch_config.json').read_text())[0]
    assert not config['provenance']['production_qm'] and not config['provenance']['exact_grafted_experiment']
    m=validate_group(mod.SCAFFOLD,config)
    assert row['inchikey']==Chem.MolToInchiKey(m)
    with pytest.raises(SystemExit) as exc:mod.main()
    assert exc.value.code==2  # preserves the review artifacts


def test_concurrent_new_motif_builds_preserve_both_manifest_merges(tmp_path):
    """Exercise real competing CLI writers, with overlap during conformer creation."""
    import csv
    import os
    import subprocess
    library=tmp_path/'library'
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    commands=[]
    for gid,frag in [('amine','[*:1]N'),('alcohol','[*:1]O')]:
        cfg=tmp_path/(gid+'.py')
        cfg.write_text("SCAFFOLD='Cc1ccc(C[*:1])cc1'\nGROUPS="+repr([
            dict(id=gid,name=gid,family='test',frag=frag,states=[('neu',0,1,0),('red',-1,2,-1)])]))
        # A deliberate overlap would expose the old read-before-build lost update.
        code="import time; from redox.build import candidates as c; original=c.build_group; " \
             "c.build_group=lambda *a: (time.sleep(.3),original(*a))[1]; c.main()"
        commands.append([sys.executable,'-c',code,'--config',str(cfg),'--library',str(library)])
    processes=[subprocess.Popen(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True) for cmd in commands]
    for p in processes:
        stdout,stderr=p.communicate(timeout=60)
        assert p.returncode==0,(stdout,stderr)
    rows=list(csv.DictReader(open(library/'manifest.csv')))
    assert len(rows)==4 and {r['id'] for r in rows}=={'amine','alcohol'}
    assert not list(library.glob('.manifest-*.tmp'))
    assert not (library/'manifest.csv.tmp').exists()


def test_geometry_identity_separately_rejects_enantiomer_and_alkene_inversion(tmp_path):
    from rdkit.Chem import AllChem
    from redox.core.common import geometry_stereo_status,smiles_to_xyz_map
    for smi,opposite in [('C[C@H](F)CC','C[C@@H](F)CC'),('F/C=C/Cl','F/C=C\\Cl')]:
        m=Chem.AddHs(Chem.MolFromSmiles(smi));assert AllChem.EmbedMolecule(m,randomSeed=43)==0
        p=tmp_path/'m.xyz';p.write_text(Chem.MolToXYZBlock(m))
        assert smiles_to_xyz_map(opposite,p) is not None  # same graph, wrong stereochemistry
        assert geometry_stereo_status(smi,p)=='ok'
        assert geometry_stereo_status(opposite,p)=='stereochemistry_mismatch'
        # Optimization outputs can carry ASE force columns; ignore them for identity.
        lines=p.read_text().splitlines()
        p.write_text('\n'.join(lines[:2]+[line+' 0 0 0' for line in lines[2:]])+'\n')
        assert geometry_stereo_status(smi,p)=='ok'


def _uma_cache(tmp_path):
    from rdkit.Chem import AllChem
    from redox.qm import uma
    smi='CCO';m=Chem.AddHs(Chem.MolFromSmiles(smi));AllChem.EmbedMolecule(m,randomSeed=9)
    seed=tmp_path/'conf_00.xyz';seed.write_text(Chem.MolToXYZBlock(m))
    xyz=tmp_path/'relaxed.xyz';xyz.write_text(seed.read_text())
    row=dict(id='m',state='neu',smiles=smi,charge='0',mult_hint='1',n_e='0')
    identity=dict(model='fake',checkpoint_sha256='abc',reference_sha256={'atom_refs':'xyz'},loader_version='1',task_name='omol')
    record=dict(cache_schema=uma.CACHE_SCHEMA,model='fake',model_identity=identity,seed_identity=uma.seed_identity([seed]),
                geometry_sha256=uma.file_sha256(xyz),energy_eV=-2.,fmax=.01,converged=True,
                n_conf=1,conf_idx=0,chosen_mult=1,mults_scanned=[1,3],n_elec=sum(a.GetAtomicNum() for a in m.GetAtoms()),
                mult_hint=1,n_e=0,id='m',state='neu',charge=0,smiles=smi)
    return uma,record,row,xyz,seed,identity


def test_uma_cache_rejects_changed_weights_seeds_geometry_hint_and_force_requirement(tmp_path):
    import copy
    uma,record,row,xyz,seed,identity=_uma_cache(tmp_path)
    assert not uma.cache_errors(record,row,xyz,[seed],identity,.02)
    assert uma.cache_errors(record,row,xyz,[seed],dict(identity,checkpoint_sha256='new'),.02)
    assert uma.cache_errors(record,dict(row,mult_hint='3'),xyz,[seed],identity,.02)
    assert uma.cache_errors(record,row,xyz,[seed],identity,.005)
    assert uma.cache_errors({k:v for k,v in record.items() if k!='cache_schema'},row,xyz,[seed],identity,.02)
    seed.write_text(seed.read_text().replace('\n\n','\nnew seed identity\n'))
    assert uma.cache_errors(record,row,xyz,[seed],identity,.02)
    record['seed_identity']=uma.seed_identity([seed])
    xyz.write_text(xyz.read_text().replace('\n\n','\nchanged output\n'))
    assert uma.cache_errors(record,row,xyz,[seed],identity,.02)


def test_uma_scan_never_selects_an_unconverged_lower_energy(tmp_path,monkeypatch):
    from ase import Atoms
    from redox.qm import uma
    monkeypatch.setattr('ase.io.read',lambda path:Atoms('H2',positions=[[0,0,0],[0,0,.7]]))
    def relax(calc,atoms,q,m,fmax,steps):
        # First seed has deceptively low energy and high residual force.
        ci=relax.calls;relax.calls+=1
        return (-100.,.5,steps,atoms) if ci==0 else (-float(m),.001,3,atoms)
    relax.calls=0;monkeypatch.setattr(uma,'relax_one',relax)
    result=uma.relax_state_scan(None,[tmp_path/'one',tmp_path/'two'],0,1,.02,5)
    assert result['chosen_mult']==3 and result['energy_eV']==-3
    assert result['converged']


def test_new_candidate_route_and_source_notes_reach_scorecard(monkeypatch):
    from redox.screening import scorecard
    group=dict(id='new',chemistry_notes='selectivity unresolved',provenance=dict(
        route_status='proposed chloride transfer',source='https://example.org/primary',exact_grafted_experiment=False))
    monkeypatch.setattr(scorecard,'_cfg',lambda name:SimpleNamespace(GROUPS=[group]))
    note=scorecard._candidate_chem_flags()['new']
    assert 'selectivity unresolved' in note and 'proposed chloride transfer' in note
    assert 'https://example.org/primary' in note and 'no verified experimental measurement' in note


def test_spin_gap_never_mixes_native_and_legacy_thermal_definitions(tmp_path,monkeypatch):
    from redox.core import common,protocol
    from redox.screening import scorecard
    from redox.validation import spins
    from redox.qm import dft
    sd=tmp_path/'spincheck/m/red_m3';sd.mkdir(parents=True)
    (sd/'opt.xyz').write_text('2\n\nH 0 0 0\nH 0 0 .7\n')
    (sd/'result.json').write_text(json.dumps(dict(charge=0,mult=3,g_thermal_eV=10,n_imag=0)))
    address=protocol.record_path(sd,0,3);address.parent.mkdir()
    address.write_text(json.dumps(dict(protocol=protocol.ACTIVE_SP,protocol_hash=address.stem,
        geom_sha1=protocol.geom_sha1(sd/'opt.xyz'),charge=0,mult=3,converged_smd=True,converged_gas=True,e_smd_eV=-9,e_gas_eV=-8)))
    monkeypatch.setattr(scorecard,'SPINCHECK',tmp_path/'spincheck');monkeypatch.setattr(scorecard,'UMA',tmp_path/'uma')
    chosen=dict(charge=0,mult=1,e_smd_eV=-10,g_thermal_eV=20,n_imag=0,
        thermal_definition='native_solution_dft_qrrho_spin',thermal_qc='ok',thermal_protocol_hash='native')
    monkeypatch.setattr(common,'read_result',lambda *a,**kw:chosen)
    monkeypatch.setattr(spins,'snapshot',lambda *args:{'geometry':'current'})
    monkeypatch.setattr(spins,'read_comparison',lambda *args:dict(status='spin_gap_verified',
        chosen_identity={'geometry':'current'},alternative_mult=3,electronic_gap_eV=1))
    monkeypatch.setattr(dft,'_opt_protocol_mismatch',lambda *args:'')
    gap,source=scorecard._min_spin_gap_eV('m')
    assert gap==1 and source.startswith('dft_electronic')  # not the mixed thermal gap of -9eV


def test_stale_uma_archive_requires_execute_and_preserves_original_bytes(tmp_path,monkeypatch):
    spec=importlib.util.spec_from_file_location('archive_test',ROOT/'scripts/pipeline/stage_uma_recalculation.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    sd=tmp_path/'uma/m/neu';sd.mkdir(parents=True)
    original='{"energy_eV": -7, "legacy": true}\n';(sd/'result.json').write_text(original)
    (sd/'relaxed.xyz').write_text('2\nold geometry\nH 0 0 0\nH 0 0 .7\n')
    monkeypatch.setattr(mod,'ROOT',tmp_path);monkeypatch.setattr(mod,'UMA',tmp_path/'uma');monkeypatch.setattr(mod,'LIBRARY',tmp_path/'library')
    monkeypatch.setattr(mod,'read_manifest',lambda:[dict(id='m',state='neu')])
    monkeypatch.setattr(mod,'model_identity',lambda model:{'model':model})
    monkeypatch.setattr(mod,'cache_errors',lambda *a:['legacy_missing_identity'])
    archive=tmp_path/'archive'
    monkeypatch.setattr(sys,'argv',['archive','--only','m','--archive-root',str(archive)])
    mod.main()
    assert sd.exists() and not archive.exists()
    monkeypatch.setattr(sys,'argv',['archive','--only','m','--archive-root',str(archive),'--execute'])
    mod.main()
    assert not sd.exists() and (archive/'m/neu/result.json').read_text()==original
    inventory=json.loads((archive/'inventory.json').read_text())
    assert inventory['status']=='archived_for_recalculation' and inventory['states'][0]['sha256']


def test_deferred_controller_preparation_never_starts_a_process(tmp_path,monkeypatch):
    spec=importlib.util.spec_from_file_location('controller_prepare_test',ROOT/'scripts/pipeline/deferred_candidate_campaign.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    monkeypatch.setattr(mod.subprocess,'Popen',lambda *a,**k:pytest.fail('preparation cannot launch any subprocess'))
    g=_group('pilot');monkeypatch.setattr(mod,'load_config',lambda name:SimpleNamespace(SCAFFOLD='Cc1ccc(C[*:1])cc1',GROUPS=[g]))
    out=tmp_path/'queue';monkeypatch.setattr(sys,'argv',['controller','--out',str(out),'--ids','pilot'])
    mod.main()
    state=json.loads((out/'status.json').read_text())
    assert state['status']=='queued_not_running' and state['quantum_jobs_submitted_at_preparation']==0
    assert state['quantum_execution_status']=='not_started'
    assert state['max_total_gpu_workers']==2 and state['threads']==2


def test_controller_waits_for_both_native_Fc_states_not_one_or_legacy(monkeypatch):
    spec=importlib.util.spec_from_file_location('controller_fc_test',ROOT/'scripts/pipeline/deferred_candidate_campaign.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    good=dict(e_smd_eV=-10,g_thermal_eV=.2,n_imag=0,thermal_definition='native_solution_dft_qrrho_spin',
              thermal_qc='ok',thermal_protocol_hash='current',sp_status='ok')
    pair={'ox':dict(good),'neu':dict(good)}
    monkeypatch.setattr(mod,'read_result',lambda gid,state:pair[state])
    assert mod.native_fc_pair_qc()[0]
    pair['neu']=dict(good,thermal_qc='missing_validated_native_thermal',g_thermal_eV=None)
    assert not mod.native_fc_pair_qc()[0]
    pair['neu']=dict(good,thermal_definition='legacy_gfn2_gas_rrho',thermal_protocol_hash=None)
    assert not mod.native_fc_pair_qc()[0]
