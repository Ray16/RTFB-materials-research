"""Durable, explicitly unqualified geometry evaluation checkpoints."""
import json
import time
from redox.core.io import atomic_write_json


def write_geometry(mol,path,comment):
    coordinates=mol.atom_coords(unit='Angstrom')
    lines=[str(mol.natm),comment]
    lines += [f'{mol.atom_pure_symbol(i)} {p[0]:.12f} {p[1]:.12f} {p[2]:.12f}' for i,p in enumerate(coordinates)]
    temp=path.with_name(path.name+'.tmp');temp.write_text('\n'.join(lines)+'\n');temp.replace(path)


def checkpoint_callback(directory):
    """SCF evaluations include rejected trials; they are not admitted optima."""
    import numpy as np
    directory.mkdir(parents=True,exist_ok=True)
    counter=0
    def callback(environment):
        nonlocal counter
        counter+=1
        mol=environment['mol'];gradient=environment['gradients']
        gradient=np.asarray(gradient.get() if hasattr(gradient,'get') else gradient)
        if not np.isfinite(gradient).all():raise ValueError('Nonfinite checkpoint gradient')
        record=dict(evaluation=counter,energy_Ha=float(environment['energy']),
                    gradient_max_Ha_Bohr=float(abs(gradient).max()),
                    SCF_converged=bool(environment['g_scanner'].converged),
                    optimizer_acceptance='not established by evaluation callback',
                    solution_geometry_qualified=False,updated=time.time())
        write_geometry(mol,directory/'last_evaluated.xyz',f'Unqualified SCF evaluation {counter}')
        with (directory/'evaluated_trajectory.xyz').open('a') as stream:
            stream.write((directory/'last_evaluated.xyz').read_text())
        with (directory/'evaluations.jsonl').open('a') as stream:
            stream.write(json.dumps(record)+'\n')
        atomic_write_json(directory/'latest.json',record)
    return callback


