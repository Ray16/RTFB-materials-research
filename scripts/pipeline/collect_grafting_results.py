#!/usr/bin/env python
"""Collect qualified reaction components and refresh the maintained visual report."""
from __future__ import annotations

import csv
import fcntl
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import time

from redox.core.common import ROOT, EV_KJ
from redox.core.io import atomic_write_json
from redox.workflows.checkpoints import verify_state_checkpoint

BASE = ROOT/'results/grafting/campaigns/2026-10-19_mecn_reference'
OUT = ROOT/'results/grafting'


def verify_sources(mapping):
    return all(Path(p).exists() and hashlib.sha256(Path(p).read_bytes()).hexdigest() == h
               for p, h in mapping.items())


def reaction_delta_g(reaction, product, components):
    """Two products minus two reactants, with explicit chemistry-specific leaving species."""
    leaving = 'chloride' if reaction['reaction_class'] == 'N_quaternization' else 'hydrogen_chloride'
    keys = (reaction['precursor_id'], 'merrifield_local_chloride', leaving)
    missing = [key for key in keys if key not in components]
    if product is None: missing.append('product_'+reaction['id'])
    if missing: return None, missing
    value = product + components[leaving] - components[reaction['precursor_id']] - components['merrifield_local_chloride']
    if not math.isfinite(value): raise ValueError('Nonfinite reaction energy')
    return value, []


def main():
    with (OUT/'.collection.lock').open('a+') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: return
        plan = json.loads((BASE/'plan.json').read_text())
        if not verify_sources(plan['input_sha256']):
            raise RuntimeError('Frozen reaction inventory or reuse evidence changed; preserve outputs for review')
        inventory = json.loads((OUT/'reactions.json').read_text())
        cached = json.loads((BASE/'reused_products.json').read_text())
        components, evidence, invalid = {}, {}, []
        for receipt in sorted((BASE/'done').glob('*.json')):
            try:
                record = verify_state_checkpoint(receipt)
                value = record['G_eV']
                if type(value) not in (int, float) or not math.isfinite(value):
                    raise ValueError('Nonfinite or missing component G')
                components[record['task']['key']] = value
                evidence[record['task']['key']] = str(receipt)
            except (ValueError, KeyError, OSError) as exc:
                invalid.append(dict(receipt=str(receipt), error=repr(exc)))
        verified_cached = {key: val for key, val in cached.items() if verify_sources(val['source_sha256'])}
        failures = {p.stem: json.loads(p.read_text()) for p in (BASE/'failures').glob('*.json')}
        rows = []
        for r in inventory['reactions']:
            product = verified_cached.get(r['id'], {}).get('G_eV', components.get('product_'+r['id']))
            value, missing = reaction_delta_g(r, product, components)
            failed_dependencies = [key for key in missing if key in failures]
            product_source = verified_cached.get(r['id'], {}).get('state_dir',
                str(Path(plan['raw_root'])/('product_'+r['id'])/'rest'))
            rows.append(dict(id=r['id'], name=r['name'], family=r['family'],
                reaction_class=r['reaction_class'], precursor_id=r['precursor_id'],
                solvent='acetonitrile', temperature_K=298.15,
                deltaG_reference_eV=value,
                deltaG_reference_kJ_mol=value*EV_KJ if value is not None else None,
                deltaG_with_base_kJ_mol=None, base_correction_required=r['base_correction_required'],
                status='qualified_reference_reaction' if value is not None else (
                    'held_failed_component_QC' if failed_dependencies else 'pending_components'),
                missing_components=';'.join(missing), failed_components=';'.join(failed_dependencies),
                atom_balanced=r['atom_balanced'], charge_balanced=r['charge_balanced'],
                product_state_path=product_source, product_smiles=r['product_smiles'],
                precursor_smiles=r['precursor_smiles']))
        table = io.StringIO(newline='')
        writer = csv.DictWriter(table, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
        tmp = OUT/'current.csv.tmp'; tmp.write_text(table.getvalue()); tmp.replace(OUT/'current.csv')
        workers = [json.loads(p.read_text()) for p in (BASE/'workers').glob('*.json')]
        claims = [p.name for p in (BASE/'claims').iterdir() if p.is_dir()]
        status = dict(updated_unix=time.time(), status='running' if claims else 'awaiting_or_terminal_workers',
            candidate_reactions=len(rows), balanced_reactions=len(rows), unique_precursors=72, shared_references=3,
            existing_product_free_energies=len(verified_cached), products_needing_QM=plan['reaction_count']-len(cached),
            completed_grafting_deltaG=sum(r['deltaG_reference_eV'] is not None for r in rows),
            qualified_components=len(components), expected_new_components=len(plan['tasks']),
            preserved_component_failures=len(failures), invalidated_receipts=invalid,
            active_claims=claims, worker_reports=len(workers),
            workers_reported_computing=sum(w.get('status')=='computing' for w in workers),
            claim_liveness='Claim presence is not proof of a live process; inspect start_verification.json',
            actual_calculation_start_records=len(list((BASE/'tasks').glob('*/status.json'))),
            new_quantum_jobs_submitted=len(plan['tasks']) if list(BASE.glob('launch_*.json')) else 0,
            queue_workers_submitted=sum(len(json.loads(p.read_text())['workers']) for p in BASE.glob('launch_*.json')),
            solvent='acetonitrile', temperature_K=298.15, base='unspecified; HCl-reference for NH/OH motifs',
            base_coupled_results_pending=True, campaign=str(BASE.relative_to(ROOT)),
            raw_root=plan['raw_root'], components=evidence,
            figure='results/grafting/figures/deltaG_ranking.pdf', browser='results/grafting/index.html',
            exact_polymer_grafting_accuracy_validated=False,
            scope='Local capped-site reference thermodynamics; no barriers/yields or base correction')
        if len(components)+len(failures)==len(plan['tasks']) and not claims:
            status['status']='complete' if not failures and not invalid else 'finished_with_QC_holds'
        atomic_write_json(OUT/'current.json', dict(status=status, candidates=rows))
        atomic_write_json(OUT/'status.json', status)
        command = [sys.executable, str(ROOT/'scripts/plotting/grafting/update_grafting_figures.py')]
        subprocess.run(command, check=True)
        print(json.dumps({k:status[k] for k in ('completed_grafting_deltaG','qualified_components','preserved_component_failures','actual_calculation_start_records')}))


if __name__=='__main__': main()
