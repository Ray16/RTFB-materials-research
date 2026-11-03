from redox import cli
from redox.workflows.catalog import catalogue
from redox.core.common import ROOT
from pathlib import Path
import pytest


def test_catalog_covers_public_commands_and_is_read_only():
    data=catalogue()
    assert {r['command'].split()[1] for r in data['commands']}==set(cli.COMMANDS)
    assert next(r for r in data['commands'] if r['command']=='redox motifs')['category']=='production'
    assert next(r for r in data['commands'] if r['command']=='redox native-thermal')['category']=='calculation_stages'
    assert any(r['kind']=='historical' for r in data['specialized'])


def test_scientific_controls_and_recovery_are_discoverable_separately():
    controls=catalogue('controls');recovery=catalogue('recovery')
    assert not controls['commands'] and not recovery['commands']
    assert all(r['category']=='controls' for r in controls['specialized'])
    assert any(r.get('script')=='scripts/fleet/run_corrected_native_queue.py' for r in controls['specialized'])
    assert any(r.get('module')=='redox.workflows.counterion_geometry_resources' for r in controls['specialized'])
    assert any(r.get('module')=='redox.workflows.charged_hessian_recovery' for r in recovery['specialized'])
    assert all(r['maintenance']=='maintained' for r in recovery['specialized'])
    assert catalogue('historical')['specialized'][0]['maintenance']=='compatibility_alias'
    with pytest.raises(ValueError):catalogue('unknown')


def test_every_discovered_advanced_command_has_an_existing_source_without_importing_it():
    commands=[]
    for row in catalogue()['specialized']:
        source=(ROOT/row['script'] if row.get('script') else ROOT/'src'/Path(row['module'].replace('.','/')).with_suffix('.py'))
        assert source.is_file(),row['command']
        commands.append(row['command'])
    assert len(commands)==len(set(commands))
