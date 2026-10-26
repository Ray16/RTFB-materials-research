#!/usr/bin/env python
"""Publish the complete independently checked frozen local-interaction panel."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--replay',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    a=parser.parse_args();replay=a.replay.resolve();out=a.out.resolve()
    if out.exists():raise FileExistsError('Preserve existing published evidence')
    d=json.loads(replay.read_text())
    if d['completed_triplets']!=12 or d['declared_triplets']!=12 or d['checked_electronic_geometries']!=36:
        raise ValueError('The panel requires twelve complete independently checked triplets')
    if any(r['missing_states'] or r['interaction_energy_kJ_mol_by_medium'] is None for r in d['states']):
        raise ValueError('A triplet is incomplete')
    if d['free_energies_computed'] or d['electronic_coupling_computed'] or d['primary_candidate_admission']:
        raise ValueError('Incorrect scientific scope')
    for path,h in d['source_sha256'].items():
        if sha(path)!=h:raise ValueError('Source changed: '+path)
    grouped={}
    for row in d['states']:
        key=row['parent_key'];grouped.setdefault(key,[]).append(row)
    if len(grouped)!=4:raise ValueError('Four distinct parent geometries required')
    summaries=[]
    for key,rows in sorted(grouped.items()):
        if sorted(round(r['cap_length_A'],2) for r in rows)!=[1.07,1.09,1.11]:
            raise ValueError('Each pose needs the three matched cap lengths')
        vals=[r['interaction_energy_kJ_mol_by_medium']['acetonitrile'] for r in rows]
        if not all(math.isfinite(v) for v in vals):raise ValueError('Nonfinite interaction')
        summaries.append(dict(parent_key=key,contact=rows[0]['parent_contact'],
            mECN_interaction_kJ_mol_min=min(vals),mECN_interaction_kJ_mol_max=max(vals),
            cap_span_kJ_mol=max(vals)-min(vals)))
    raw_hash=sha(replay)
    report=dict(status='all_twelve_frozen_local_interaction_triplets_independently_reconstructed',
        declared_triplets=12,completed_triplets=12,checked_electronic_geometries=36,
        states=d['states'],per_parent_cap_sensitivity=summaries,
        definition=d['definition'],scope=d['scope'],source_replay=str(replay),
        source_replay_sha256=raw_hash,source_sha256=d['source_sha256'],
        publisher_sha256=sha(__file__),original_collector_status=d['status'],
        status_explanation='The earlier collector uses a fixed pending status string; all input states and all twelve computed sums were checked here.',
        basis_superposition_correction_computed=False,free_energies_computed=False,
        electronic_coupling_computed=False,conformer_convergence_established=False,
        primary_candidate_admission=False)
    out.mkdir(parents=True);(out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'completed':12,'summaries':summaries},indent=2))


if __name__=='__main__':main()
