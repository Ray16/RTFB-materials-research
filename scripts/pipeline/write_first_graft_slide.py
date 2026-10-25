#!/usr/bin/env python
"""Generate a scoped result slide from independently reconstructed first-site G."""
import argparse
import json
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve(); report = json.loads(source.read_text())
    api.verify_sources(report['source_sha256'])
    if report['status'] != 'available_first_site_grafting_independently_reconstructed' or report['synthesis_admission']:
        raise ValueError('Scoped, independently reconstructed first-site evidence is required')
    rows = []
    for row in sorted(report['singles'],key=lambda r:r['site']):
        values = row['by_cutoff']['50']
        rows.append(f"{row['site']} & ${values['first_graft_kJ_mol']['free_ions']:.1f}$ & "
            f"${values['first_graft_kJ_mol']['contact_pair']:.1f}$ & ${values['polymer_minus_capped_kJ_mol']:+.2f}$" + r'\\')
    text = '\n'.join([
        r'\begin{frame}{PMDI first-site attachment: a sampled polymer model}',
        r'\small Five-repeat open chain; SMD MeCN, 298.15 K; DBU-assisted reactions.',
        r'\medskip\centering',
        r'\begin{tabular}{lrrr}\toprule',
        r'Site & Separated ions & Contact pair & Polymer minus local\\',
        r' & \multicolumn{3}{c}{kJ mol$^{-1}$}\\\midrule',
        *rows, r'\bottomrule\end{tabular}',
        r'\par\medskip\raggedright',
        r'\textbf{Context penalty:} common precursor and base terms cancel from the polymer-minus-local comparison.',
        r'\par\medskip\textbf{Scope:} one sampled minimum per occupied graph. Site preference changes sign across the 25--100 cm$^{-1}$ rotor settings.',
        r'\par\medskip\textbf{Pending:} paired-state free energy, converged conformer ensembles, selectivity and accessible utilization.',
        r'\par\smallskip\footnotesize Independent state/nuclear reconstruction, 25 October 2026. No material candidate is admitted.',
        r'\end{frame}', ''])
    out = ROOT/'slides/generated/pmdi_first_graft_2026-10-25.tex';out.write_text(text)
    atomic_write_json(out.with_suffix('.json'),dict(source=str(source),source_sha256=api.digest(source),
        producer_sha256=api.digest(__file__),tex_sha256=api.digest(out),synthesis_admission=False,paired_native_G_available=False))


if __name__=='__main__':main()
