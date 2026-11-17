#!/usr/bin/env python3
"""Focused slide figure from the frozen first corrected occupied batch analysis."""
import json
import os
from pathlib import Path


def main():
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from redox.core.common import ROOT
    from redox.core.hashing import sha256_file
    folder = ROOT/'slides/evidence/2026-11-02_occupied_sampling'
    source = folder/'finite_batches.json'; data = json.loads(source.read_text())
    cell = next(c for c in data['cells'] if c['model_key'] == 'all_same__occ_2_3')
    selected = cell['partitions']['0.1']; cutoffs = ['25', '50', '100']
    values = [selected[k]['comparison']['ensemble_G_change_kJ_mol'] for k in cutoffs]
    limit = selected['50']['comparison']['tolerance_kJ_mol']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 18, 'axes.labelsize': 19,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, ax = plt.subplots(figsize=(7.6, 4.9))
    ax.axvspan(-limit, 0, color='#EDF4F3', zorder=0)
    ax.barh([2, 1, 0], values, height=.50, color=['#537E83', '#007F86', '#537E83'])
    ax.axvline(-limit, color='#BD6928', linestyle='--', linewidth=2)
    ax.text(-limit-.05, 2.59, f'−{limit:.2f} tolerance', ha='right', color='#96501D', fontsize=17)
    for y, value in zip([2, 1, 0], values):
        ax.text(value-.07, y, f'−{abs(value):.2f}', ha='right', va='center', fontsize=20, color='#15343B')
    ax.set(yticks=[2, 1, 0], yticklabels=cutoffs, xlim=(-3.5, 0), ylim=(-.65, 2.95),
           xlabel=r'$G_{20}-G_{10}$ (kJ/mol)', ylabel=r'Rotor cutoff (cm$^{-1}$)')
    ax.set_xticks([-3, -2, -1, 0]); ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout(pad=.9)
    target = ROOT/'results/figures/slides/occupied_sampling_convergence.pdf'
    target.parent.mkdir(parents=True, exist_ok=True); fig.savefig(target); plt.close(fig)
    record = dict(figure=str(target.relative_to(ROOT)), figure_sha256=sha256_file(target),
        source_sha256=sha256_file(source), producer_sha256=sha256_file(__file__),
        model_key=cell['model_key'], rotor_cutoffs_cm=cutoffs, ensemble_G_change_kJ_mol=values,
        tolerance_kJ_mol=limit, statistical_intervals_plotted=False, experimental_populations=False)
    (folder/'figure.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record))


if __name__ == '__main__': main()
