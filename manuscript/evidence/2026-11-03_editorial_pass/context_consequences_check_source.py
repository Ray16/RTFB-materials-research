"""Check the quantitative interpretation tests against the current frozen inputs."""
import json
from pathlib import Path
from redox.core.hashing import sha256_file
from redox.validation.context_consequences import assess

NAMES=('base_assisted.json','pmdi_first_site_grafting.json','pmdi_capacity.json','spacing_conformers.json','tempo_solution_raw_result.json')


def audit(folder):
    expected=json.loads((folder/'context_consequences.json').read_text())
    for name in NAMES:
        if sha256_file(folder/name)!=expected['source_files_sha256'][name]:
            raise ValueError('Context comparison input changed')
    actual=assess(*(json.loads((folder/name).read_text()) for name in NAMES))
    for name,value in actual.items():
        if expected.get(name)!=value:raise ValueError('Context interpretation does not reconstruct: '+name)
    return dict(reaction_sign_screen_changes=actual['route_definition']['changes_at_declared_nonpositive_screen'],
        complete_reactions=actual['route_definition']['complete'],
        unchanged_within_base_pairwise_orders=actual['route_definition']['unchanged_pairwise_orders'],
        positive_polymer_context_penalties=actual['polymer_context']['positive_penalties'],
        site_preference_changes_sign=actual['polymer_context']['site_preference_changes_sign'],
        experimental_or_algorithmic_superiority_inferred=False)
