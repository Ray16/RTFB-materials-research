"""Shared duplicate-aware comparison of complete nested native input sets."""
from redox.properties.conformer_ensemble import batch_change, cluster_minima, summarize


def finite_batches(rows, matrix, first_count, *, temperature_K=298.15,
                   rmsd_cutoffs_A=(.05, .1, .2), energy_tolerance_kJ_mol=.2,
                   rotor_cutoffs_cm=(25, 50, 100), tolerance_kJ_mol=2):
    """Restrict full-set well assignments to each nested batch; hits are not weights."""
    early = [row for row in rows if int(row['key'].split('conf')[-1]) <= first_count]
    early_keys = {row['key'] for row in early}
    scenarios = {}
    for threshold in rmsd_cutoffs_A:
        groups = cluster_minima(rows, matrix, rmsd_cutoff_A=threshold,
                                energy_tolerance_kJ_mol=energy_tolerance_kJ_mol)
        first_groups = [[key for key in group if key in early_keys] for group in groups]
        first_groups = [group for group in first_groups if group]
        values = {}
        for cutoff in rotor_cutoffs_cm:
            rotor = str(cutoff)
            records = [dict(key=r['key'], E_eV=r['E_eV'], G_eV=r['G_by_cutoff_eV'][rotor]) for r in rows]
            first = summarize([r for r in records if r['key'] in early_keys], first_groups,
                              [r['key'] for r in early], temperature_K)
            expanded = summarize(records, groups, [r['key'] for r in rows], temperature_K)
            values[rotor] = dict(first=first, expanded=expanded,
                                comparison=batch_change(first, expanded, tolerance_kJ_mol))
        scenarios[str(threshold)] = values
    return scenarios
