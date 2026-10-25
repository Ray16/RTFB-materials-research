from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def test_only_complete_evaluation_frames_are_retained_without_qualification(tmp_path):
    module = api.load_module('interrupted_geometry_frames_test', ROOT / 'scripts/pipeline/recover_geometry_interruptions.py')
    complete = '2\nUnqualified SCF evaluation 1\nH 0 0 0\nH 0 0 .74\n'
    path = tmp_path / 'trajectory.xyz'
    path.write_text(complete + '2\nUnqualified SCF evaluation 2\nH 0 0 0\n')
    assert module.evaluated_frames(path) == [complete]
