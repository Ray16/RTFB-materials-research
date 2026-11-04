"""One checked sampling result must survive flat evidence-prefix relocation."""
from copy import deepcopy
import json

import pytest

from redox.core.hashing import sha256_file
from redox.validation import sampling_publication as api


def fixture(tmp_path, monkeypatch, prefix):
    monkeypatch.setattr(api, 'ROOT', tmp_path)
    folder = tmp_path/'evidence'; folder.mkdir()
    def path(name): return folder/(prefix+name)
    comparisons = {k: dict(within_batch_tolerance=False, tolerance_kJ_mol=2.,
        minimum_E_change_kJ_mol=-3.88, minimum_G_change_kJ_mol=-.65,
        ensemble_G_change_kJ_mol=value) for k, value in [('25', -2.55), ('50', -2.40), ('100', -2.54)]}
    partitions = {r: {k: dict(comparison=deepcopy(c), expanded=dict(unique_minima=20))
                         for k, c in comparisons.items()} for r in ('0.05', '0.1', '0.2')}
    cell = dict(model_key='all_same__occ_2_3', checked_inputs=20, partitions=partitions)
    saved = dict(complete_cells=1, occupied_ensemble_cycles_qualified=False, cells=[cell], source_sha256={})
    path('finite_batches.json').write_text(json.dumps(saved))
    for name in ('native_snapshot.json', 'receipt_audit.json', 'sampling_contract.json'):
        path(name).write_text('{}')
    def reconstruct(*args):
        assert args == tuple(path(n) for n in ('native_snapshot.json','receipt_audit.json','sampling_contract.json'))
        return deepcopy(saved)
    monkeypatch.setattr(api, 'analyze', reconstruct)
    plotter=tmp_path/'scripts/plotting/plot_occupied_sampling.py';plotter.parent.mkdir(parents=True);plotter.write_text('fixture plot producer')
    figure=tmp_path/'figure.pdf';figure.write_bytes(b'fixture figure')
    path('figure.json').write_text(json.dumps(dict(figure='figure.pdf',figure_sha256=sha256_file(figure),
        producer_sha256=sha256_file(plotter),source_sha256=sha256_file(path('finite_batches.json')),
        ensemble_G_change_kJ_mol=[-2.55,-2.40,-2.54],tolerance_kJ_mol=2.,
        statistical_intervals_plotted=False,experimental_populations=False)))
    macros={'SamplingTwentyFiveAbs':'2.55','SamplingFiftyAbs':'2.40','SamplingHundredAbs':'2.54',
            'SamplingTolerance':'2.00','SamplingMinimumEDrop':'3.88'}
    path('values.tex').write_text(''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k,v in macros.items()))
    path('table.tex').write_text(api.sampling_table(cell))
    manifest=dict(status='frozen_first_corrected_occupied_sampling_cell',model_key=cell['model_key'],
        files={n:sha256_file(path(n)) for n in ('finite_batches.json','native_snapshot.json',
               'receipt_audit.json','sampling_contract.json','figure.json','values.tex')})
    path('manifest.json').write_text(json.dumps(manifest))
    return folder, path, manifest


@pytest.mark.parametrize('prefix', ['', 'occupied_sampling_'])
def test_deck_and_flat_manuscript_epoch_use_the_same_checked_quantities(tmp_path, monkeypatch, prefix):
    folder, _, _ = fixture(tmp_path, monkeypatch, prefix)
    result=api.audit(folder,prefix)
    assert result['ensemble_G_change_kJ_mol']['50']==-2.40
    assert result['distinct_minima']==20
    assert not result['within_batch_tolerance'] and not result['occupied_ensemble_cycles_qualified']


@pytest.mark.parametrize('bad', ['table', 'macro', 'plot', 'scope'])
def test_changed_numbers_or_promoted_scope_cannot_reach_publication(tmp_path, monkeypatch, bad):
    folder,path,manifest=fixture(tmp_path,monkeypatch,'occupied_sampling_')
    if bad=='table':path('table.tex').write_text(path('table.tex').read_text().replace('-3.88','-0.38'))
    if bad=='macro':path('values.tex').write_text(path('values.tex').read_text().replace('2.40','0.24'))
    if bad=='plot':
        data=json.loads(path('figure.json').read_text());data['statistical_intervals_plotted']=True
        path('figure.json').write_text(json.dumps(data))
    if bad=='scope':
        data=json.loads(path('finite_batches.json').read_text());data['occupied_ensemble_cycles_qualified']=True
        path('finite_batches.json').write_text(json.dumps(data))
    # A publisher updating file hashes still cannot change the numerical meaning.
    manifest['files']={n:sha256_file(path(n)) for n in manifest['files']}
    path('manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError):api.audit(folder,'occupied_sampling_')
