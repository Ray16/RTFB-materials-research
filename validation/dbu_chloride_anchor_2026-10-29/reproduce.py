#!/usr/bin/env python3
"""Conditional conversion of a published association constant; no fitting or QM."""
import hashlib
import json
import math
from pathlib import Path

ARCHIVE = Path('/nfs/lambda_stor_01/homes/rzhu/archive/0_redox')
SOURCE = ARCHIVE / 'results/grafting/base_assisted/current.json'
EXPECTED_SHA256 = 'd26186a8e587e849c426b0ca40b9761907ce71e317a7d3ffabff7d4a79872d20'
raw = SOURCE.read_bytes()
if hashlib.sha256(raw).hexdigest() != EXPECTED_SHA256:
    raise SystemExit('Archived computed source changed; review before reuse')
data = json.loads(raw)
correction = data['corrections']['dbu']
model = correction['ion_pair_association_kJ_mol']
assert abs((correction['contact_pair_kJ_mol'] - correction['free_ions_kJ_mol']) - model) < 1e-9
R = 8.31446261815324
T = 298.15  # comparison assumption, not verified as the experimental temperature
K_M_inv = 419.0  # reported concentration association constant; no uncertainty supplied here
reference_M = 1.0
experimental_conditional = -R * T * math.log(K_M_inv * reference_M) / 1000
model_K_M_inv = math.exp(-model * 1000 / (R * T)) / reference_M
print(json.dumps(dict(
    status='preliminary_conditional_comparison',
    source_path=str(SOURCE), source_sha256=EXPECTED_SHA256,
    temperature_K_assumed=T, standard_concentration_mol_L=reference_M,
    experimental_reported_K_M_inv=K_M_inv,
    experimental_derived_deltaG_kJ_mol=experimental_conditional,
    model_association_deltaG_kJ_mol=model, model_implied_K_M_inv=model_K_M_inv,
    model_minus_experiment_kJ_mol=model-experimental_conditional,
    limitations=['experimental temperature and error not verified',
                 'CD3CN versus modeled MeCN, PF6 and ionic activities not matched',
                 'experiment apparent binding versus modeled one contact-pair state',
                 'not a full grafting reaction validation or transferable correction'],
    quantum_calculations_run=False), indent=2))
