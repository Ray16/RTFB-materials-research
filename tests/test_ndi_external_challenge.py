import csv
import importlib.util
from pathlib import Path

import pytest


source = Path(__file__).resolve().parents[1] / 'scripts/validation/prepare_ndi_external_challenge.py'
spec = importlib.util.spec_from_file_location('ndi_external_challenge', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def inputs(tmp_path, doi=module.DOI, peaks='−1040 mV and −1440 mV, −880 mV and −1350 mV'):
    xml = tmp_path / 'source.xml'
    xml.write_text(f'''<article><article-id pub-id-type="doi">{doi}</article-id>
      <p>The targeted styrene–NDI monomer was prepared from 4-vinylaniline and 1,4,5,8-naphthalenetetracarboxylic dianhydride.</p>
      <p>In acetonitrile with TBAClO4 and ferrocene/ferrocenium, quasi-reversible peaks: {peaks}.</p>
      <p>calculated C30H18N2O4 (470.1267)</p></article>''')
    library = tmp_path / 'library.csv'
    with library.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['id', 'smiles']); writer.writeheader()
    search = tmp_path / 'search.csv'
    with search.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['id', 'inchikey']); writer.writeheader()
    return xml, library, search


def test_source_peaks_produce_proxies_without_benchmark_admission(tmp_path):
    result = module.prepare(*inputs(tmp_path))
    assert result['molecule']['formula'] == 'C30H18N2O4'
    assert result['molecule']['graph_checks'] == {'N_para_vinylphenyl': 2, 'imide_N': 2, 'terminal_vinyl': 2}
    assert result['measurement']['peaks'][0]['peak_midpoint_V'] == pytest.approx(-0.960)
    assert result['measurement']['peaks'][1]['peak_midpoint_V'] == pytest.approx(-1.395)
    assert result['measurement']['peaks'][0]['peak_separation_V'] == pytest.approx(0.160)
    assert result['strict_benchmark_admission'] is False and result['computed_potentials'] is None


@pytest.mark.parametrize('doi,peaks', [('10.1234/wrong', '−1040 mV and −1440 mV, −880 mV and −1350 mV'),
                                       (module.DOI, '−1040 mV and −1440 mV, −880 mV')])
def test_source_mismatch_is_rejected(tmp_path, doi, peaks):
    with pytest.raises(ValueError):
        module.prepare(*inputs(tmp_path, doi=doi, peaks=peaks))
