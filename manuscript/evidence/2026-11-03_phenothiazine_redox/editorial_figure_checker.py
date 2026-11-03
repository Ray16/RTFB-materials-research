"""Verify that the final figure layout preserves every frozen numerical value."""
import json
from redox.core.hashing import sha256_file

def audit(folder,figures):
    proof=json.loads((folder/'editorial_figure_review.json').read_text())
    source=folder/'motif_definition_comparison.json'
    if sha256_file(source)!=proof['source_sha256']:raise ValueError('Editorial source differs')
    data=json.loads(source.read_text());rows=[next(r for r in data['states'] if r['parent_key']==k) for k in proof['keys']]
    if ([r['whole_chain_nonadditivity_E_kJ_mol'] for r in rows]!=proof['whole_chain_E_kJ_mol']
            or [r['frozen_methyl_motif_interaction_E_kJ_mol'] for r in rows]!=proof['motif_E_kJ_mol']
            or proof['error_bars_or_population_weights_added']):raise ValueError('Editorial numeric interpretation changed')
    for name in ('viologen_interaction_definitions_final.pdf','viologen_interaction_definitions_final_grayscale.pdf'):
        expected=next(h for p,h in proof['output_sha256'].items() if p.endswith('/'+name))
        if sha256_file(figures/name)!=expected:raise ValueError('Final figure differs from its checked render')
    if sha256_file(folder/'editorial_figure_producer.py')!=proof['producer_sha256']:raise ValueError('Figure producer changed')
    return dict(frames=4,checked_values=8,numerical_values_unchanged=True,new_calculations=False,
                color_and_grayscale_figures_bound=True)
