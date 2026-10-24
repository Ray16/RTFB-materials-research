"""The public default must actually reach current engines and keep shard ownership."""
from pathlib import Path

import pytest

from redox import current
from redox.workflows import launch_current


@pytest.mark.parametrize('module',['uma','dft'])
def test_current_commands_preserve_shards_cpus_and_GPU_gate(module,tmp_path):
    jobs=launch_current.commands(module,[1,3],tmp_path,threads=2,shard_count=5,shard_start=2,
        cpu_sets=[[4,5],[6,7]],extra=['--only','design_pmdi_ethyl'])
    assert len(jobs)==2
    for index,job in enumerate(jobs):
        assert job[job.index('--gpu')+1]==str((1,3)[index])
        assert job[job.index('--cpus')+1]==('4,5','6,7')[index]
        assert job[job.index('--shard')+1]==f'5:{index+2}'
        assert 'redox.workflows.gpu' in job and 'redox.qm.uma_current' in job
        assert job[job.index('--stage')+1]==module
        assert job[-2:]==['--only','design_pmdi_ethyl']


def test_public_commands_are_updated_only_during_dispatch(monkeypatch):
    original=current.original.COMMANDS
    def inspect(argv):
        assert current.original.COMMANDS['uma'][0]=='redox.qm.uma_current'
        assert current.original.COMMANDS['launch'][0]=='redox.workflows.launch_current'
        assert current.original.COMMANDS['dft'][1]==('--stage','dft')
        return 17
    monkeypatch.setattr(current.original,'main',inspect)
    assert current.main(['uma'])==17
    assert current.original.COMMANDS is original
