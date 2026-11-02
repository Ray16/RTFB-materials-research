"""Electronic-only implementation contributions to declared balanced cycles."""
import math
from redox.core.common import HARTREE_EV, EV_KJ


def reconstruct(package, energies):
    states={r['id']:r for r in package['states']}
    if len(states)!=len(package['states']) or set(energies)!=set(states):
        raise ValueError('Every exact package state must be checked once')
    differences={}
    for key,row in energies.items():
        values=[row[k] for k in ('stored_Ha','installed_Ha','corrected_Ha')]
        if any(type(v) not in (int,float) or not math.isfinite(v) for v in values):
            raise ValueError('Finite electronic energies required')
        if (row['charge'],row['mult'])!=(states[key]['charge'],states[key]['mult']):
            raise ValueError('Electronic state differs from declared package')
        old,installed,corrected=values
        differences[key]=dict(installed_minus_stored_Ha=installed-old,
            corrected_minus_stored_Ha=corrected-old,corrected_minus_installed_Ha=corrected-installed)
    rows=[]
    for reaction in package['reactions']:
        coeff=reaction['coefficients']
        if not coeff or not set(coeff)<=set(states) or any(type(c) not in (int,float) or not math.isfinite(c) for c in coeff.values()):
            raise ValueError('Finite coefficients of explicit exact states required')
        if sum(c*states[k]['charge'] for k,c in coeff.items())!=0:
            raise ValueError('Actual cycle must conserve total charge')
        # Atom balance is established from the frozen package's numerical replay;
        # this scalar sensitivity does not reassert that physical qualification.
        if reaction['kind']=='one_electron_potential_vs_Fc':
            if reaction['electrons']!=1:raise ValueError('Declared one-electron referenced cycle required')
            factor=-HARTREE_EV;unit='V'
        elif reaction['kind'] in ('grafting_free_energy','reaction_free_energy'):
            factor=HARTREE_EV*EV_KJ;unit='kJ/mol'
        else:raise ValueError('Unsupported reaction quantity')
        shifts={name.replace('_Ha',''):math.fsum(coeff[k]*differences[k][name] for k in coeff)*factor
                for name in next(iter(differences.values()))}
        rows.append(dict(id=reaction['id'],kind=reaction['kind'],unit=unit,
                         electronic_contribution_to_shift=shifts))
    return dict(status='declared_cycle_electronic_sensitivity_reconstructed',states=len(states),cycles=rows,
        nuclear_terms_recomputed=False,total_free_energy_requalified=False,
        new_potentials_or_reaction_G_qualified=False)
