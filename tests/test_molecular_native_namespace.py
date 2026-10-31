import csv
import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from redox.core.protocol import ACTIVE_SP,geom_sha1,protocol_hash,record_path
from redox.workflows.molecular_native_namespace import source_state,safe_component,command


@pytest.fixture
def archive(tmp_path):
    root=tmp_path/'archive';directory=root/'calcs/dft/molecule/neu';directory.mkdir(parents=True)
    mol=Chem.AddHs(Chem.MolFromSmiles('CCO'));AllChem.EmbedMolecule(mol,randomSeed=17)
    Chem.MolToXYZFile(mol,str(directory/'opt.xyz'))
    raw=dict(id='molecule',state='neu',charge=0,mult=1,optimized=True,opt_xc='r2scan',opt_basis='def2-svp',opt_disp='d4',solvent='acetonitrile')
    (directory/'result.json').write_text(json.dumps(raw))
    (root/'library').mkdir()
    with (root/'library/manifest.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=['id','state','charge','mult_hint','smiles']);writer.writeheader()
        writer.writerow(dict(id='molecule',state='neu',charge=0,mult_hint=1,smiles='CCO'))
    path=record_path(directory,0,1);path.parent.mkdir()
    record=dict(protocol=ACTIVE_SP,protocol_hash=path.stem,geom_sha1=geom_sha1(directory/'opt.xyz'),charge=0,mult=1,
        e_smd_eV=-100.,e_gas_eV=-99.,converged_smd=True,converged_gas=True,
        source='computed:unit_fixture',software={'pyscf':'2.14.0','gpu4pyscf':'1.8.1'},smd_spin_contam=0.)
    path.write_text(json.dumps(record));return root,directory,path


def test_exact_source_geometry_and_matched_active_energy_are_admitted(archive):
    root,directory,path=archive
    row,raw,energy=source_state(root,'molecule','neu','CCO')
    assert energy==path and raw['mult']==1 and row['smiles']=='CCO'


@pytest.mark.parametrize('change',['charge','geometry_protocol','active_basis','spin_nonfinite','hydrogen','partial_native'])
def test_incompatible_source_or_partial_native_outputs_are_preserved(archive,change):
    root,directory,path=archive
    raw=json.loads((directory/'result.json').read_text());record=json.loads(path.read_text())
    if change=='charge':raw['charge']=1
    if change=='geometry_protocol':raw['opt_basis']='def2-svpd'
    if change=='active_basis':record['protocol']=dict(ACTIVE_SP,basis='def2-tzvp')
    if change=='spin_nonfinite':record['smd_spin_contam']=float('nan')
    if change=='hydrogen':
        lines=(directory/'opt.xyz').read_text().splitlines();lines[0]=str(int(lines[0])-1)
        (directory/'opt.xyz').write_text('\n'.join(lines[:-1])+'\n')
    if change=='partial_native':(directory/'native_thermal').mkdir()
    (directory/'result.json').write_text(json.dumps(raw));path.write_text(json.dumps(record))
    with pytest.raises((ValueError,FileExistsError)):source_state(root,'molecule','neu','CCO')


@pytest.mark.parametrize('identifier',['','../escape','/absolute','nested/name','.','..'])
def test_identifiers_cannot_escape_the_private_state_tree(identifier):
    with pytest.raises(ValueError):safe_component(identifier)


def test_only_private_dft_locks_repairs_and_case_paths_are_writable(tmp_path):
    case=tmp_path/'case';sources=tmp_path/'sources';case.mkdir();sources.mkdir()
    spec=dict(prior_root='/original',archive='/archive',case=str(case),sources=str(sources))
    args=command(spec)
    binds=[tuple(args[i+1:i+3]) for i,v in enumerate(args) if v=='--bind']
    assert (str(case/'dft'),'/original/calcs/dft') in binds
    assert (str(case/'locks'),'/original/logs/pipeline_locks') in binds
    assert not any(target in ['/archive','/original'] for _,target in binds)
