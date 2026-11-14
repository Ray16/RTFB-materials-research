#!/usr/bin/env python3
"""Show minimum and ensemble changes in two separately sampled backbones."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def plot(source, output):
    data = json.loads(source.read_text())
    if (data['status'] != 'source_bound_finite_native_batch_analysis'
            or data['global_ensemble_convergence_proven'] or data['configurations_are_pooled']):
        raise ValueError('Separate finite checked backbone analyses required')
    if any(sha256_file(path) != digest for path, digest in data['source_sha256'].items()):
        raise ValueError('Batch analysis input or producer changed')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
    figure, axis = plt.subplots(figsize=(7.8, 4.4), layout='constrained')
    metrics = ['minimum_E_change_kJ_mol', 'minimum_G_change_kJ_mol', 'ensemble_G_change_kJ_mol']
    records = []
    for i, (group, color, marker) in enumerate(zip(data['groups'], ['#087F8C', '#C76E00'], ['o', '^'])):
        if group['status'] != 'finite_native_batches_compared' or group['occupied_sites']:
            raise ValueError('Complete empty-backbone input sets required')
        nominal = group['partitions']['0.1']['50']['comparison']
        x = [j+(-.11 if i == 0 else .11) for j in range(3)]
        values, ranges = [], []
        for j, metric in enumerate(metrics):
            settings = [rotor['comparison'][metric] for rmsd in group['partitions'].values()
                        for rotor in rmsd.values()]
            values.append(nominal[metric])
            ranges.append([min(settings), max(settings)])
            axis.plot([x[j], x[j]], ranges[-1], color=color, lw=2.5, alpha=.65, zorder=2)
            axis.plot([x[j]-.025, x[j]+.025], [ranges[-1][0]]*2, color=color, lw=1)
            axis.plot([x[j]-.025, x[j]+.025], [ranges[-1][1]]*2, color=color, lw=1)
        axis.scatter(x, values, color=color, marker=marker, s=62,
                     label=group['configuration'].replace('_', ' '), zorder=3)
        records.append(dict(configuration=group['configuration'], nominal=nominal,
                            settings_ranges_kJ_mol=dict(zip(metrics, ranges)),
                            unique_minima=group['partitions']['0.1']['50']['expanded']['unique_minima']))
    axis.axhline(0, color='#34454C', lw=.8)
    axis.axhline(-2, color='#65727A', ls='--', lw=1)
    axis.text(2.35, -1.96, '2 kJ/mol batch criterion', ha='right', va='bottom', fontsize=10, color='#53626B')
    axis.set(xticks=range(3), xticklabels=['Lowest E', 'Lowest G', 'Ensemble G'],
             xlim=(-.45, 2.45), ylim=(-2.25, .30),
             ylabel='Change from 10 to 20 inputs (kJ/mol)')
    axis.grid(axis='y', alpha=.12)
    axis.legend(loc='lower left', frameon=False, ncol=2)
    output.parent.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for extension in ('pdf', 'png'):
        path = output.with_suffix('.'+extension)
        figure.savefig(path, dpi=240, facecolor='white')
        outputs[str(path.resolve())] = sha256_file(path)
    plt.close(figure)
    record = dict(status='finite_empty_backbone_batch_stability_figure', records=records,
                  analysis=str(source.resolve()), analysis_sha256=sha256_file(source),
                  source_sha256=sha256_file(Path(__file__)), outputs_sha256=outputs,
                  nominal_rotor_cm=50, nominal_duplicate_RMSD_A=.1,
                  settings_ranges_are_confidence_intervals=False,
                  global_ensemble_convergence_proven=False,
                  caption='Points: 50 cm^-1 rotor and 0.1 A duplicate cutoff. Ranges: explicit rotor/duplicate settings; not uncertainty intervals. Each backbone configuration is a separate species.')
    atomic_write_json(output.with_suffix('.json'), record)
    return record


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    plot(args.analysis.resolve(), args.out.resolve())
