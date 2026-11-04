"""An interaction label must identify its actual audited parent/cap inputs."""
import json
from pathlib import Path
import pytest
from redox.core.hashing import sha256_file
from redox.workflows.motif_interaction_energy import validate_construction_binding


def setup(tmp_path):
    tasks=[]
    for label,charge,atoms in [('A',2,28),('B',2,28),('AB',4,56)]:
        xyz=tmp_path/(label+'.xyz');graph=tmp_path/(label+'.mol');xyz.write_text('explicit input '+label);graph.write_text('explicit graph '+label)
        tasks.append(dict(id=label,geometry=str(xyz),molecular_graph=str(graph),geometry_sha256=sha256_file(xyz),
            graph_sha256=sha256_file(graph),charge=charge,smiles='test graph '+label,atoms=atoms))
    entry=dict(key='parent__cap109',parent_key='parent',cap_length_A=1.09,molecules=tasks)
    (tmp_path/'construction.json').write_text(json.dumps(dict(entries=[entry])))
    return dict(sources=str(tmp_path)),dict(source_id=entry['key'],parent_key='parent',cap_length_A=1.09,tasks=tasks),dict(parent_key='parent',cap_length_A=1.09)


def test_correct_parent_and_copied_input_bindings_are_retained(tmp_path):
    validate_construction_binding(*setup(tmp_path))


@pytest.mark.parametrize('change',['parent','cap','geometry','graph'])
def test_labels_or_inputs_cannot_be_swapped_after_construction(tmp_path,change):
    spec,plan,report=setup(tmp_path)
    if change=='parent':report['parent_key']='another_parent'
    if change=='cap':report['cap_length_A']=1.11
    if change=='geometry':Path(plan['tasks'][0]['geometry']).write_text('different frozen coordinates')
    if change=='graph':Path(plan['tasks'][1]['molecular_graph']).write_text('different graph')
    with pytest.raises(ValueError):validate_construction_binding(spec,plan,report)
