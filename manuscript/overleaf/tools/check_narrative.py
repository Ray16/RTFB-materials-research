#!/usr/bin/env python3
"""Check the fixed benchmark and the surfaced matched-parent comparisons."""
from pathlib import Path
import collections
import csv
import hashlib
import json
import re

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA256 = '1a7f35323aee304e2eeb017ea8a24ed0fd79cd65c31d8a10d0c15ff8eaa908bf'
SELECTED = [
    ('tempo_parent', '1', r'TEMPO, first oxidation\cite{gerken2015tempo}'),
    ('bm_n_ethylphenothiazine', '1', r'N-Ethylphenothiazine, first oxidation\cite{attanayake2019}'),
    ('bm_1_benzyl_1_methyl_4_4_bipyridinium', '1', r'Benzyl--methyl viologen, first reduction\cite{kim2022}'),
    ('bm_1_benzyl_1_methyl_4_4_bipyridinium', '2', r'Benzyl--methyl viologen, second reduction\cite{kim2022}'),
]

def audit():
    source = ROOT / 'data/benchmark_endpoints.csv'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == SOURCE_SHA256, 'Benchmark source changed'
    rows = list(csv.DictReader(source.open()))
    counts = dict(collections.Counter(r['family'] for r in rows))
    assert counts == {'phenothiazine': 5, 'quinone': 22, 'nitroxide': 1, 'aromatic imide': 9, 'viologen': 2}
    lookup = {(r['id'], r['wave']): r for r in rows}
    assert len(lookup) == len(rows) == 39
    errors = [float(r['calculated_V_vs_Fc']) - float(r['experimental_V_vs_Fc']) for r in rows]
    assert all(abs(e - float(r['signed_error_V'])) < 1e-12 for e, r in zip(errors, rows))
    means = {wave: sum(float(r['signed_error_V']) for r in rows if r['wave'] == wave) / sum(r['wave'] == wave for r in rows) for wave in ['1', '2']}
    assert sum(r['wave'] == '1' for r in rows) == 30
    assert all(float(r['signed_error_V']) > 0 for r in rows if r['wave'] == '1')
    assert f"{means['1']:+.3f}" == '+0.266' and f"{means['2']:+.3f}" == '-0.064'
    paired = []
    for r in rows:
        if r['wave'] == '2':
            first = lookup[(r['id'], '1')]
            calculated = abs(float(r['calculated_V_vs_Fc']) - float(first['calculated_V_vs_Fc']))
            measured = abs(float(r['experimental_V_vs_Fc']) - float(first['experimental_V_vs_Fc']))
            assert calculated > measured
            paired.append({'id': r['id'], 'calculated_spacing_V': calculated, 'experimental_spacing_V': measured})
    assert len(paired) == 9
    text = (ROOT / 'main.tex').read_text()
    selected = []
    for ident, wave, label in SELECTED:
        r = lookup[(ident, wave)]
        values = [float(r[k]) for k in ['experimental_V_vs_Fc', 'calculated_V_vs_Fc', 'signed_error_V']]
        expected = label + ' & ' + ' & '.join(f'${v:+.3f}$' for v in values) + r'\\'
        assert expected in text, 'Matched-parent table does not match source: ' + ident + '/' + wave
        selected.append({'id': ident, 'wave': wave, 'experimental_V': values[0], 'calculated_V': values[1], 'signed_error_V': values[2]})
    original = json.loads((ROOT / 'collaboration_origin.json').read_text())
    for name, digest in original['initial_export_members'].items():
        if name.startswith('inputs/'):
            assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, 'Frozen numerical input changed: ' + name
    combined = text + (ROOT / 'supporting_information.tex').read_text()
    definitions = re.findall(r'\\(?:newcommand|renewcommand)\{\\([A-Za-z]+)\}', '\n'.join(x.read_text() for x in (ROOT / 'inputs').glob('*.tex')))
    for document in ['main.tex', 'supporting_information.tex']:
        content = (ROOT / document).read_text()
        labels = re.findall(r'\\label\{([^}]+)\}', content)
        assert len(labels) == len(set(labels)), 'Duplicate labels: ' + document
        assert set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', content)) <= set(labels), 'Undefined labels: ' + document
        for name in re.findall(r'\\(?:includegraphics)(?:\[[^\]]*\])?\{([^}]+)\}', content):
            assert (ROOT / name).is_file(), 'Missing figure: ' + name
    citations = set()
    for group in re.findall(r'\\cite(?:[a-zA-Z]*)?(?:\[[^\]]*\])?\{([^}]+)\}', combined):
        citations.update(group.split(','))
    bib = (ROOT / 'references.bib').read_text() + (ROOT / 'novelty_references.bib').read_text()
    assert citations <= set(re.findall(r'@\w+\{([^,]+),', bib)), 'Missing bibliography keys'
    assert 'Claude Code' in text and 'OpenAI Codex' in text
    assert 'tree/refactor/clean-workflow' not in text
    return {'status': 'passed', 'benchmark_source_sha256': SOURCE_SHA256, 'endpoints': len(rows), 'family_counts': counts, 'MAE_V': sum(abs(e) for e in errors)/len(errors), 'bias_V': sum(errors)/len(errors), 'mean_signed_error_by_step_V': means, 'paired_spacings': paired, 'matched_parent_rows': selected, 'frozen_numerical_inputs_unchanged': True, 'references_and_figure_paths_checked': True}

if __name__ == '__main__':
    result = audit()
    print(json.dumps(result, indent=2))
