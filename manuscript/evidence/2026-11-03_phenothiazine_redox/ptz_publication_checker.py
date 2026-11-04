"""Replay the frozen finite phenothiazine oxidation result for paper and slides."""
import json
import math

from redox.core.hashing import sha256_file
from redox.properties.reference_redox import neutral_oxidation_vs_fc


def values(report):
    cycles=report['cycles']; primary=cycles['minus']['potential_by_cutoff_V']['50']
    all_values=[v for row in cycles.values() for v in row['potential_by_cutoff_V'].values()]
    difference=abs(primary-cycles['plus']['potential_by_cutoff_V']['50'])*1e6
    return (f'\\newcommand{{\\PTZPrimaryV}}{{{primary:.3f}}}\n'
            f'\\newcommand{{\\PTZLowV}}{{{min(all_values):.3f}}}\n'
            f'\\newcommand{{\\PTZHighV}}{{{max(all_values):.3f}}}\n'
            f'\\newcommand{{\\PTZBranchMicroV}}{{{difference:.2f}}}\n')


def audit(folder):
    report=json.loads((folder/'ptz_redox_cycle.json').read_text())
    if (set(report['cycles'])!={'minus','plus'} or not report['redox_potential_qualified']
            or report['conformer_ensemble_qualified'] or report['experimental_accuracy_validated']
            or report['material_admission']):
        raise ValueError('Finite qualified oxidation scope required')
    for branch,cycle in report['cycles'].items():
        if neutral_oxidation_vs_fc(cycle['components'])!=cycle:
            raise ValueError('Frozen phenothiazine cycle differs from replay')
        for role,component in cycle['components'].items():
            leaf=branch if role=='graft_neu' else 'oxidized' if role=='graft_ox' else role
            path=folder/('ptz_review_'+leaf+'.json'); review=json.loads(path.read_text())
            if (sha256_file(path) not in component['sources'].values()
                    or review['G_by_cutoff_eV']!=component['G_by_cutoff_eV']
                    or (review['charge'],review['mult'])!=(component['state']['charge'],component['state']['multiplicity'])
                    or not math.isfinite(review['minimum_frequency_cm']) or review['minimum_frequency_cm'] < -5
                    or not 0<=review['maximum_column_error']<=.01
                    or review['actual_SCF_checkpoint_frames_checked']!=component['actual_native_SCF_frames_checked']):
                raise ValueError('Frozen component review does not bind its value or QC')
    if (folder/'ptz_redox_values.tex').read_text()!=values(report):
        raise ValueError('Phenothiazine TeX values differ')
    return dict(status='frozen_phenothiazine_redox_arithmetic_and_component_scope_checked',
                primary_V=report['cycles']['minus']['potential_by_cutoff_V']['50'],
                branches=2, rotor_cutoffs=3, new_quantum_validation=False,
                experimental_accuracy_validated=False, material_admission=False)
