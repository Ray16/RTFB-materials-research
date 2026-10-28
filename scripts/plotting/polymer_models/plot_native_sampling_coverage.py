#!/usr/bin/env python3
"""Plot source-checked conformer energies without imputing missing minima."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.properties.conformer_ensemble import EV_KJ_MOL


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--replay', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    data = json.loads(a.replay.read_text())
    if (data['status'] != 'archived_polymer_native_physical_records_replayed'
            or data['quantum_calculations_run'] or data['global_ensemble_convergence_proven']):
        raise ValueError('Checked finite native replay required')
    for path, digest in data['producer_sha256'].items():
        if sha256_file(path) != digest:
            raise ValueError('Replay producer changed')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.8), sharex=True, sharey=True,
                             layout='constrained')
    summaries = []
    for i, configuration in enumerate(('all_same', 'alternating')):
        for j, sites in enumerate(([], [2])):
            ax = axes[i, j]
            rows = [r for r in data['records'] if r['configuration'] == configuration
                    and r['occupied_sites'] == sites]
            if not rows or len({r['graph_sha256'] for r in rows}) != 1:
                raise ValueError('One exact graph per panel required')
            rows.sort(key=lambda r: int(r['key'].split('conf')[-1]))
            indices = [int(r['key'].split('conf')[-1]) for r in rows]
            if len(set(indices)) != len(indices) or not set(indices) <= set(range(1, 21)):
                raise ValueError('Distinct declared indices 1–20 required')
            values = [r['G_by_cutoff_eV']['50'] for r in rows]
            minimum = min(values)
            relative = [(value-minimum)*EV_KJ_MOL for value in values]
            missing = sorted(set(range(1, 21))-set(indices))
            ax.scatter(indices, relative, color='#087F8C', s=43, label='Checked native G', zorder=3)
            ax.scatter(missing, [-.9]*len(missing), marker='x', color='#6B747B', s=38,
                       label='Missing native G', clip_on=False)
            ax.axvline(10.5, color='#BFC7CA', ls=':', lw=1)
            ax.axhline(0, color='#34454C', lw=.7)
            name = 'Empty backbone' if not sites else 'Site 2 occupied'
            ax.set_title(f'{configuration.replace("_", " ")} · {name} ({len(rows)}/20)', loc='left', fontsize=12)
            ax.set_xlim(.5, 20.5)
            ax.set_ylim(-1.8, 15)
            ax.set_xticks([1, 5, 10, 15, 20])
            ax.grid(axis='y', alpha=.12)
            if i == 1:
                ax.set_xlabel('Declared conformer input index')
            if j == 0:
                ax.set_ylabel('G − lowest checked G (kJ/mol)')
            summary = dict(configuration=configuration, occupied_sites=sites,
                           checked=len(rows), missing_indices=missing,
                           observed_range_kJ_mol=max(relative),
                           lowest_checked_key=rows[values.index(minimum)]['key'])
            earlier = [value for index, value in zip(indices, values) if index <= 10]
            if earlier:
                summary['lowest_checked_minus_lowest_checked_in_first_ten_kJ_mol'] = (minimum-min(earlier))*EV_KJ_MOL
            summaries.append(summary)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside lower center', ncol=2, frameon=False)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for extension in ('pdf', 'png'):
        destination = a.out.with_suffix('.'+extension)
        fig.savefig(destination, dpi=240, facecolor='white')
        outputs[str(destination.resolve())] = sha256_file(destination)
    plt.close(fig)
    atomic_write_json(a.out.with_suffix('.json'), dict(
        status='finite_checked_native_sampling_coverage', panels=summaries,
        source_replay=str(a.replay.resolve()), replay_sha256=sha256_file(a.replay),
        producer_sha256=sha256_file(Path(__file__)), outputs_sha256=outputs,
        model='alternate smooth-surface acetonitrile', temperature_K=298.15,
        rotor_cutoff_cm=50, global_convergence_proven=False,
        caption=('Each panel compares one exact graph/configuration only. Missing records are '
                 'marked below the zero line, not assigned physical energies. Incomplete batches '
                 'cannot establish a 10-to-20 convergence comparison. No polymer reaction or '
                 'experimental occupancy probability is inferred.')))


if __name__ == '__main__':
    main()
