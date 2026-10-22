#!/usr/bin/env python
"""Plot repeated real descriptor timings without relabeling them as quantum speed."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--out', type=Path, default=ROOT / 'results/figures/validation')
    args = parser.parse_args()
    source = json.loads(args.summary.read_text())
    if source['status'] != 'real_descriptor_concurrency_checked' or source['serial_discrepancies']:
        raise ValueError('Only independently checked descriptor runs may be plotted')
    rows = source['conditions']
    if any(row['repeats'] != 3 for row in rows):
        raise ValueError('This figure declares three repeats per condition')
    sources = {str(p):sha(p) for p in (args.summary, Path(__file__))}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,
                         'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 3.7))
    figure.subplots_adjust(left=.09, right=.98, top=.80, bottom=.28, wspace=.43)
    for hosts, color, label in [(1,'#0072B2','One host'),(2,'#D55E00','Two hosts')]:
        selected = sorted([row for row in rows if len(row['hosts']) == hosts], key=lambda r:r['workers'])
        if not selected:
            continue
        x = np.array([row['workers'] for row in selected], dtype=float)
        if hosts == 2:
            x += .22
        for axis, key in zip(axes, ('launch_to_observed_completion_seconds','workers_executing_chunks')):
            values = [row[key]['median'] for row in selected]
            errors = [[row[key]['median']-row[key]['minimum'] for row in selected],
                      [row[key]['maximum']-row[key]['median'] for row in selected]]
            axis.errorbar(x, values, yerr=errors, fmt='o-' if hosts == 1 else 's', color=color,
                          capsize=3, lw=1, markersize=4, label=label)
    for axis, letter in zip(axes, ('a','b')):
        axis.set(xticks=[1,2,4,8],xlabel='Requested CPU workers',xlim=(.5,8.7))
        axis.text(-.17,1.06,letter,transform=axis.transAxes,fontweight='bold',fontsize=10)
    axes[0].set(ylabel='Launch to observed completion (s)',ylim=(0,None))
    axes[1].set(ylabel='Workers executing at least one chunk',ylim=(0,9),yticks=[0,2,4,6,8])
    axes[1].plot([1,8],[1,8],color='.65',ls='--',lw=.7,zorder=0)
    axes[0].legend(frameon=False,fontsize=7)
    figure.suptitle('Measured CPU descriptor concurrency: 837 exact designs',y=.98,fontsize=11)
    figure.text(.09,.87,'Same 1,593 recipes · 8 tasks per chunk · 3 repeats per condition',fontsize=8)
    figure.text(.09,.16,'Elapsed time includes node preflight, SSH, scorer setup and work; preparation is excluded.',fontsize=7)
    figure.text(.09,.105,'Ranges are observed min–max, not confidence intervals. Shared nodes and a small molecular library.',fontsize=7)
    figure.text(.09,.05,'This measures descriptor execution; it does not measure million-molecule chemistry or quantum throughput.',fontsize=7)
    args.out.mkdir(parents=True, exist_ok=True)
    stem = args.out / 'cpu_descriptor_concurrency'
    outputs = []
    for suffix in ('pdf','png','svg'):
        destination = stem.with_suffix('.'+suffix)
        figure.savefig(destination,dpi=600);outputs.append(destination)
    plt.close(figure)
    flat = []
    for row in rows:
        flat.append(dict(condition=row['condition'],workers=row['workers'],hosts=len(row['hosts']),repeats=row['repeats'],
                         **{key+'_'+stat:value for key in ('launch_to_observed_completion_seconds','worker_timestamp_span_seconds','workers_executing_chunks','largest_worker_peak_RSS_KiB') for stat,value in row[key].items()}))
    with stem.with_suffix('.csv').open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(flat[0]));writer.writeheader();writer.writerows(flat)
    outputs.append(stem.with_suffix('.csv'))
    if any(sha(path) != expected for path,expected in sources.items()):
        raise ValueError('Figure input changed during plotting')
    atomic_write_json(stem.with_suffix('.json'),dict(status='checked_descriptor_concurrency_figure',
                      source_sha256=sources,output_sha256={str(p):sha(p) for p in outputs},
                      quantum_throughput=False,experimental_accuracy=False))


if __name__ == '__main__':
    main()
