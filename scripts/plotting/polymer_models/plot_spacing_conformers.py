#!/usr/bin/env python
"""Render accepted finite-control conformers with their actual contact geometry."""
import argparse
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from rdkit import Chem

from redox.build.polymer_cases import _path
from redox.build.polymer_contacts import contact_metrics, pi_core_groups
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

COLORS = ('#D55E00', '#CC79A7')
BACKBONE = '#0072B2'
EV_KJ = 1.602176634e-19 * 6.02214076e23 / 1000  # Exact SI definitions.


def load_checked_rows(path):
    review = json.loads(path.read_text())
    if review['checked_native_G'] != 12 or review['pending_keys']:
        raise ValueError('All twelve independently reconstructed finite-baseline states are required')
    api.verify_sources(review['source_sha256'])
    renderer_path = ROOT/'scripts/plotting/polymer_models/oligomer_summary.py'
    renderer = api.load_module('exact_saved_molecular_mapping', renderer_path)
    states = {row['key']: row for row in review['states']}
    rows = []
    sources = {str(p.resolve()):api.digest(p) for p in (path, Path(__file__), renderer_path,
        ROOT/'src/redox/build/polymer_contacts.py', ROOT/'src/redox/build/polymer_cases.py')}
    for graph in ('occ_2_3', 'occ_2_4'):
        pair = [states[f'{graph}__pose{i}'] for i in (0, 1)]
        minimum = min(s['G_eV'] for s in pair)
        for state in pair:
            matches = [Path(p) for p in review['source_sha256'] if '/thermal/tasks/'+state['key']+'/status.json' in p]
            if len(matches) != 1:
                raise ValueError('Unique native state report required')
            report_path = matches[0]
            report = json.loads(report_path.read_text())
            task = report['task']
            xyz_path = Path(state['state_dir'])/'opt.xyz'
            mol, mapping = renderer.mapped_molecule(dict(task, final_xyz=str(xyz_path)))
            mapped_groups = [[mapping[i] for i in group] for group in pi_core_groups(mol, task['fragment'])]
            if {frozenset(g) for g in mapped_groups} != {frozenset(g) for g in task['pi_groups']}:
                raise ValueError('The conjugated cores do not map to the declared XYZ atoms')
            mapped = Chem.MolFromSmiles(task['mapped_smiles'])
            canonical = Chem.MolFromSmiles(task['smiles'])
            match = canonical.GetSubstructMatch(mapped, useChirality=True)
            if len(match) != mapped.GetNumAtoms():
                raise ValueError('Mapped backbone and canonical stereochemical graph differ')
            backbone = [match[i] for _, i in sorted(_path(mapped, task['repeat_units']).items())]
            coords = np.asarray(mol.GetConformer().GetPositions())
            # Coordinates are now in canonical graph order, independently of
            # the XYZ atom order. Contact metrics use that same graph order.
            groups = pi_core_groups(mol, task['fragment'])
            metric = contact_metrics(coords, groups)[0]
            from ase.io import read
            original_metric = contact_metrics(read(xyz_path).positions, task['pi_groups'])[0]
            for key in ('centroid_distance_A','interplane_angle_deg','plane_separation_A','nearest_core_atom_A'):
                if abs(metric[key]-original_metric[key]) > 1e-8:
                    raise ValueError('Mapped drawing changed a contact metric')
            heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
            center = coords[heavy].mean(axis=0)
            _, _, rotation = np.linalg.svd(coords[heavy]-center, full_matrices=False)
            if np.linalg.det(rotation) < 0:
                rotation[-1] *= -1  # Preserve handedness; never reflect a conformer.
            view = (coords-center)@rotation.T
            tilt = np.deg2rad(25)
            camera = np.array([[1,0,0],[0,np.cos(tilt),-np.sin(tilt)],[0,np.sin(tilt),np.cos(tilt)]])
            view = view@camera
            view -= (view[heavy].max(axis=0)+view[heavy].min(axis=0))/2
            rows.append(dict(key=state['key'], graph=graph, task=task, mol=mol, groups=groups,
                backbone=backbone, coords=view, heavy=heavy, contact=metric,
                core_coordinates_A=[coords[g].tolist() for g in groups],
                coordinates_A=coords.tolist(), atomic_numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                coordinate_rotation=(rotation.T@camera).tolist(),
                G_eV=state['G_eV'], deltaG_conformer_kJ_mol=(state['G_eV']-minimum)*EV_KJ,
                xyz=str(xyz_path), atoms=mol.GetNumAtoms(), source_report=str(report_path)))
            sources[str(xyz_path)] = api.digest(xyz_path)
            sources[str(report_path)] = api.digest(report_path)
    return rows, sources, None


def load_cached_rows(path):
    """Render archived graph-ordered coordinates without opening cluster parents."""
    cached = json.loads(path.read_text())
    rows = cached['states']
    expected = {f'occ_{sites}__pose{pose}' for sites in ('2_3','2_4') for pose in (0,1)}
    if len(rows) != 4 or {r['key'] for r in rows} != expected:
        raise ValueError('The cached finite-control figure requires all four states')
    for row in rows:
        mol = Chem.AddHs(Chem.MolFromSmiles(row['task']['smiles']))
        if [a.GetAtomicNum() for a in mol.GetAtoms()] != row['atomic_numbers']:
            raise ValueError('Cached atom order differs from the molecular graph')
        coords = np.asarray(row['coordinates_A'], dtype=float)
        rotation = np.asarray(row['coordinate_rotation'], dtype=float)
        if coords.shape != (mol.GetNumAtoms(),3) or not np.isfinite(coords).all():
            raise ValueError('Invalid cached coordinates')
        if (rotation.shape != (3,3) or not np.allclose(rotation@rotation.T,np.eye(3),atol=1e-10,rtol=0)
                or abs(np.linalg.det(rotation)-1) > 1e-10):
            raise ValueError('The view must preserve handedness and spatial scale')
        groups = pi_core_groups(mol,row['task']['fragment'])
        if {frozenset(g) for g in groups} != {frozenset(g) for g in row['groups']}:
            raise ValueError('Cached conjugated cores differ from the molecular graph')
        metric = contact_metrics(coords, groups)[0]
        if any(abs(metric[k]-row['contact'][k]) > 1e-8 for k in ('centroid_distance_A','interplane_angle_deg','plane_separation_A')):
            raise ValueError('Cached contact metrics do not reconstruct')
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
        view = (coords-coords[heavy].mean(axis=0))@rotation
        view -= (view[heavy].max(axis=0)+view[heavy].min(axis=0))/2
        row.update(mol=mol, coords=view, heavy=heavy)
    sources = {str(p.resolve()):api.digest(p) for p in (path, Path(__file__), ROOT/'src/redox/build/polymer_contacts.py')}
    return rows, sources, cached['source_sha256']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--review', type=Path)
    source.add_argument('--cached', type=Path, help='Render saved coordinate data without original calculation reads')
    parser.add_argument('--out', type=Path, default=ROOT/'results/figures/polymer_models/spacing_conformers')
    args = parser.parse_args()
    if args.cached and any(args.out.with_suffix('.'+ext).exists() for ext in ('json','pdf','svg','png')):
        parser.error('Cached replay requires a fresh --out; preserve the original figure and its evidence')
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    rows, sources, original_sources = load_cached_rows(args.cached) if args.cached else load_checked_rows(args.review)
    # Close-core first in each row; the order is a geometric display choice,
    # not an imposed free-energy or global-minimum ranking.
    rows.sort(key=lambda r: (r['graph'], r['contact']['centroid_distance_A']))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    figure = plt.figure(figsize=(7.2, 4.4))
    horizontal_span = max(np.ptp(r['coords'][r['heavy'],0]) for r in rows)*1.12
    vertical_span = max(np.ptp(r['coords'][r['heavy'],1]) for r in rows)*1.3
    for index, row in enumerate(rows):
        axis = figure.add_subplot(2, 2, index+1)
        coords, mol = row['coords'], row['mol']
        core = {atom:COLORS[k] for k, group in enumerate(row['groups']) for atom in group}
        backbone = set(row['backbone'])
        for bond in sorted(mol.GetBonds(), key=lambda b:coords[[b.GetBeginAtomIdx(),b.GetEndAtomIdx()],2].mean()):
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if mol.GetAtomWithIdx(a).GetAtomicNum() == 1 or mol.GetAtomWithIdx(b).GetAtomicNum() == 1:
                continue
            color, width = '#A7A7A7', .85
            if a in backbone and b in backbone:
                color, width = BACKBONE, 2.5
            elif a in core and b in core and core[a] == core[b]:
                color, width = core[a], 1.9
            axis.plot(*coords[[a,b],:2].T, color=color, linewidth=width, zorder=2)
        axis.scatter(*coords[row['heavy'],:2].T,
            c=[core.get(i, BACKBONE if i in backbone else '#A7A7A7') for i in row['heavy']],
            s=[11 if i in core or i in backbone else 6 for i in row['heavy']], zorder=3)
        centroids = np.array([coords[g].mean(axis=0) for g in row['groups']])
        axis.plot(*centroids[:,:2].T, color='#303030', linestyle='--', linewidth=1, zorder=4)
        for label_index, repeat in enumerate(row['task']['occupied_sites']):
            point = coords[row['backbone'][2*repeat-1]]
            axis.annotate(str(repeat), xy=point[:2], xytext=(-12,5) if label_index == 0 else (8,-7),
                textcoords='offset points', color=BACKBONE, fontsize=8, fontweight='bold', zorder=5,
                arrowprops=dict(arrowstyle='-',color=BACKBONE,lw=.6),
                bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=.3))
        axis.set_xlim(-horizontal_span/2, horizontal_span/2)
        axis.set_ylim(-vertical_span/2, vertical_span/2)
        axis.set_aspect('equal')
        axis.set_axis_off()
        occupancy = 'Adjacent (2,3)' if row['graph'] == 'occ_2_3' else 'One spacer (2,4)'
        descriptor = 'close cores' if index % 2 == 0 else 'separated cores'
        axis.set_title(f'{chr(97+index)}  {occupancy}: {descriptor}', loc='left', fontsize=9, fontweight='bold', pad=0)
        metric = row['contact']
        axis.text(.01, -.13,
            f"Centroids {metric['centroid_distance_A']:.2f} Å; planes {metric['interplane_angle_deg']:.1f}°\n"
            f"Conformer ΔG = +{row['deltaG_conformer_kJ_mol']:.1f} kJ mol⁻¹",
            transform=axis.transAxes, fontsize=7.5, va='top')
    handles = [Line2D([0],[0],color=BACKBONE,lw=2.5,label='Open backbone'),
        Line2D([0],[0],color=COLORS[0],lw=2,label='Core 1'),
        Line2D([0],[0],color=COLORS[1],lw=2,label='Core 2'),
        Line2D([0],[0],color='#303030',linestyle='--',label='Centroid distance')]
    figure.legend(handles=handles, loc='upper center', ncol=4, frameon=False, fontsize=8)
    figure.subplots_adjust(left=.025,right=.985,bottom=.15,top=.84,wspace=.06,hspace=.85)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    for ext in ('pdf','svg','png'):
        figure.savefig(args.out.with_suffix('.'+ext),dpi=600,facecolor='white')
    plt.close(figure)
    exported = [{k:v for k,v in r.items() if k not in ('mol','coords','heavy')} for r in rows]
    caption = ('Actual SMD(MeCN) phthalimide five-repeat geometries from the twelve-state finite baseline. '
        'Adjacent sites (2,3) and one-spacer sites (2,4) each have a close-core and a separated-core local minimum. '
        'Numbers identify the occupied repeats on the open backbone; ungrafted site 3 in the one-spacer model retains benzyl chloride. '
        'Hydrogens are hidden; these orthographic projections use proper rotations of the saved XYZ coordinates and a common spatial scale. Covalent bonds follow the exact checked graph. '
        'Conformer ΔG is relative to the lower G in the same two-pose occupancy graph at 298.15 K and a 50 cm^-1 cutoff. '
        'These within-graph differences include backbone and solvent response; four-state site nonadditivity is a distinct quantity. '
        'The sampled pairs do not establish a converged ensemble, isolated core binding or an electron-transfer rate.')
    api.verify_sources(sources)
    atomic_write_json(args.out.with_suffix('.json'),dict(status='cached_coordinates_rendered_no_physical_requalification' if args.cached else 'accepted_saved_conformers_mapped_and_rendered',
        source_sha256=sources, states=exported, common_projected_horizontal_span_A=horizontal_span,
        common_projected_vertical_span_A=vertical_span, proper_rotation_only=True, caption=caption,
        cached_replay=bool(args.cached), original_source_sha256=original_sources,
        conformer_convergence_established=False, electronic_coupling_computed=False,
        output_sha256={str(args.out.with_suffix('.'+ext)):api.digest(args.out.with_suffix('.'+ext)) for ext in ('pdf','svg','png')}))
    print(json.dumps(dict(rendered=len(rows),figure=str(args.out.with_suffix('.pdf')))))


if __name__ == '__main__':
    main()
