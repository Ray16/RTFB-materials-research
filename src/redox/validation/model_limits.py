"""Scope observed model limits without changing calculations or censoring residuals."""
import hashlib
import json
from pathlib import Path

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.protocol import ACTIVE_SP


def selection_admission(smiles,charges,role,*,review_path=None):
    """Qualify decision use, separately from solver/native-record admission.

    A declared +2/+1 oxidation remains in scope even if the model's window walk
    stopped before that step. This prevents a failed second-wave model from
    turning an uncertain capacity into a confidently predicted one-electron path.
    Experimental benchmark residuals must not call this decision gate.
    """
    path=Path(review_path) if review_path is not None else ROOT/'config/model_validation_limits.json'
    try:
        review=json.loads(path.read_text())
        if review['schema']!=1 or not isinstance(review['rules'],list):raise ValueError('schema')
    except (OSError,ValueError,TypeError,KeyError) as exc:
        raise RuntimeError(f'Model-selection review unavailable: {path}') from exc
    molecule=Chem.MolFromSmiles(smiles) if isinstance(smiles,str) else smiles
    if molecule is None:raise ValueError('Model-selection scope requires a valid molecular graph')
    result=dict(model_primary_admissible=True,model_validation_status='no_scoped_limitation',
                model_validation_evidence=None,model_max_observed_error_V=None)
    for rule in review['rules']:
        query=Chem.MolFromSmiles(rule['core_smiles'])
        if query is None:raise ValueError('Invalid model-review core')
        if (role!=rule['role'] or not {rule['oxidized_charge'],rule['reduced_charge']}<=set(charges)
                or not molecule.HasSubstructMatch(query)):continue
        evidence=ROOT/rule['evidence_path']
        current=(evidence.is_file() and hashlib.sha256(evidence.read_bytes()).hexdigest()==rule['evidence_sha256'])
        if not current:status='model_review_evidence_stale'
        elif rule['tested_energy_protocol']!=ACTIVE_SP:status='model_protocol_revalidation_required'
        elif rule['status']=='limitation_observed':status=rule['id']+'_unqualified'
        else:status='unrecognized_model_review_status'
        result.update(model_primary_admissible=False,model_validation_status=status,
                      model_validation_evidence=rule['evidence_path'],
                      model_max_observed_error_V=rule['observed_potential_error_V'])
        break
    return result
