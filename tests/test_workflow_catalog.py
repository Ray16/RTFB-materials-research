from redox import cli
from redox.workflows.catalog import catalogue


def test_catalog_covers_public_commands_and_is_read_only():
    data=catalogue()
    assert {r['command'].split()[1] for r in data['commands']}==set(cli.COMMANDS)
    assert next(r for r in data['commands'] if r['command']=='redox motifs')['category']=='production'
    assert next(r for r in data['commands'] if r['command']=='redox native-thermal')['category']=='calculation_stages'
    assert any(r['kind']=='historical' for r in data['specialized'])
