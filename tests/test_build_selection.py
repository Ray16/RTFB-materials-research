"""A malformed batch selection must never fall back to building the registry."""
import sys

import pytest

from redox.build import candidates


@pytest.fixture
def definitions(tmp_path,monkeypatch):
    config=tmp_path/'candidates.py'
    config.write_text('SCAFFOLD="C[*:1]"\nGROUPS='+repr([
        dict(id='alcohol',name='alcohol',family='test',frag='[*:1]O',states=[('neu',0,1,0),('ox',1,2,1)]),
        dict(id='amine',name='amine',family='test',frag='[*:1]N',states=[('neu',0,1,0),('ox',1,2,1)])])+'\n')
    library=tmp_path/'library'
    monkeypatch.setattr(candidates,'LIBRARY',library)
    monkeypatch.setattr(candidates,'MANIFEST',library/'manifest.csv')
    return config,library


def test_batch_dry_run_selects_both_exact_ids_without_writing(definitions,monkeypatch,capsys):
    config,library=definitions
    monkeypatch.setattr(sys,'argv',['builder','--config',str(config),'--library',str(library),
                                    '--only','alcohol,amine','--dry-run'])
    candidates.main();output=capsys.readouterr().out
    assert '[valid] alcohol:' in output and '[valid] amine:' in output
    assert not library.exists()


@pytest.mark.parametrize('selection',['alcohol,missing',''])
def test_unknown_or_empty_batch_cannot_create_any_outputs(definitions,monkeypatch,selection):
    config,library=definitions
    monkeypatch.setattr(sys,'argv',['builder','--config',str(config),'--library',str(library),'--only',selection])
    with pytest.raises(SystemExit,match='invalid --only'):candidates.main()
    assert not library.exists()
