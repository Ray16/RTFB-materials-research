"""Shared source-bound occupied-sampling checks for the manuscript and slides."""
import json

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.validation.sampling_batches import analyze


def audit(folder, prefix=""):
    def artifact(name): return folder/(prefix+name)
    manifest = json.loads((artifact('manifest.json')).read_text())
    if (manifest['status'] != 'frozen_first_corrected_occupied_sampling_cell'
            or any(sha256_file(artifact(name)) != h for name, h in manifest['files'].items())):
        raise ValueError('Frozen sampling publication sources changed')
    saved = json.loads((artifact('finite_batches.json')).read_text())
    if any(sha256_file(p) != h for p, h in saved['source_sha256'].items()):
        raise ValueError('The original finite-batch input or producer changed')
    # Repeats graph/RMSD/partition analysis of checked records, not quantum QC.
    repeated = analyze(artifact('native_snapshot.json'), artifact('receipt_audit.json'), artifact('sampling_contract.json'))
    if {k: v for k, v in repeated.items() if k != 'source_sha256'} != {
            k: v for k, v in saved.items() if k != 'source_sha256'}:
        raise ValueError('Finite-cell geometry or partition reconstruction differs')
    if saved['complete_cells'] != 1 or saved['occupied_ensemble_cycles_qualified']:
        raise ValueError('Only the first complete cell, without balanced-cycle admission, is shown')
    cell = next(c for c in saved['cells'] if c['model_key'] == manifest['model_key'])
    primary = cell['partitions']['0.1']['50']; comparison = primary['comparison']
    if any(r['comparison']['within_batch_tolerance'] for cutoff in cell['partitions'].values() for r in cutoff.values()):
        raise ValueError('The first occupied sampling publication requires its failed batch diagnostic')
    changes = {k: cell['partitions']['0.1'][k]['comparison']['ensemble_G_change_kJ_mol'] for k in ('25', '50', '100')}
    macros = {'SamplingTwentyFiveAbs': f"{abs(changes['25']):.2f}",
              'SamplingFiftyAbs': f"{abs(changes['50']):.2f}",
              'SamplingHundredAbs': f"{abs(changes['100']):.2f}",
              'SamplingTolerance': f"{comparison['tolerance_kJ_mol']:.2f}",
              'SamplingMinimumEDrop': f"{abs(comparison['minimum_E_change_kJ_mol']):.2f}"}
    figure = json.loads((artifact('figure.json')).read_text())
    root = ROOT
    if (figure['figure_sha256'] != sha256_file(root/figure['figure'])
            or figure['source_sha256'] != sha256_file(artifact('finite_batches.json'))
            or figure['producer_sha256'] != sha256_file(root/'scripts/plotting/plot_occupied_sampling.py')
            or figure['ensemble_G_change_kJ_mol'] != [changes[k] for k in ('25', '50', '100')]
            or figure['tolerance_kJ_mol'] != comparison['tolerance_kJ_mol']
            or figure['statistical_intervals_plotted'] or figure['experimental_populations']):
        raise ValueError('Sampling chart data or its interpretation differs')
    expected = ''.join('\\newcommand{\\'+k+'}{'+v+'}\n' for k, v in macros.items())
    if (artifact('values.tex')).read_text() != expected:
        raise ValueError('Published numbers differ from the checked finite comparison')
    if artifact('table.tex').exists() and artifact('table.tex').read_text() != sampling_table(cell):
        raise ValueError('Sampling table differs from the checked nested comparison')
    return dict(model_key=cell['model_key'], checked_inputs=cell['checked_inputs'],
        distinct_minima=primary['expanded']['unique_minima'],
        ensemble_G_change_kJ_mol=changes, minimum_E_change_kJ_mol=comparison['minimum_E_change_kJ_mol'],
        tolerance_kJ_mol=comparison['tolerance_kJ_mol'], within_batch_tolerance=False,
        RMSD_and_partition_analysis_repeated=True, original_physical_QC_repeated=False,
        occupied_ensemble_cycles_qualified=False, numbers=list(macros.values()))


def sampling_table(cell):
    """Exact rounded SI values; no uncertainty or convergence inferred."""
    rows = [r"\begin{tabular}{rrrr}\toprule",
            r"Rotor cutoff / cm$^{-1}$ & $\Delta E_{\min}$ & $\Delta G_{\min}$ & $\Delta G_{\mathrm{set}}$\\\midrule"]
    for cutoff in ('25', '50', '100'):
        c = cell['partitions']['0.1'][cutoff]['comparison']
        rows.append(f"{cutoff} & {c['minimum_E_change_kJ_mol']:.2f} & {c['minimum_G_change_kJ_mol']:.2f} & {c['ensemble_G_change_kJ_mol']:.2f}" + r"\\")
    rows.append(r"\bottomrule\end{tabular}")
    return '\n'.join(rows)+'\n'
