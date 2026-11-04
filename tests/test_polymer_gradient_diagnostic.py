import json

import pytest

from redox.workflows.polymer_gradient_diagnostic import select_points


def test_only_failed_column_base_and_signed_displacements_are_selected(tmp_path):
    checks=[dict(atom=0,axis=1,step_Bohr=.0005,relative_error=.02)]
    (tmp_path/'derivative_checks.json').write_text(json.dumps(dict(checks=checks)))
    for i,x in enumerate([0.,-.0005,.0005]):
        folder=tmp_path/'SCF_evaluations'/f'{i:04d}';folder.mkdir(parents=True)
        (folder/'evaluation.json').write_text(json.dumps(dict(coordinates_Bohr=[[0,x,0]])))
    failed,paths=select_points(tmp_path)
    assert failed==checks[0] and [p.parent.name for p in paths]==['0000','0001','0002']
    checks[0]['relative_error']=.001;(tmp_path/'derivative_checks.json').write_text(json.dumps(dict(checks=checks)))
    with pytest.raises(ValueError):select_points(tmp_path)
