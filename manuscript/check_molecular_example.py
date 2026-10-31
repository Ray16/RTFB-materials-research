"""Check the frozen exact-molecule CPU example against the versioned raw data."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def audit(folder,root=None):
    path=folder/'exact_molecular_CPU_validation.json'
    if not path.exists():return None
    root=Path(root) if root is not None else Path(__file__).resolve().parents[1]
    validation=json.loads(path.read_text());frozen=json.loads((folder/'exact_molecular_CPU_replay.json').read_text())
    package=root/'examples/exact_viologen_replay/data';manifest=json.loads((package/'manifest.json').read_text())
    if (validation['status']!='extracted_native_cache_fresh_numpy_replay_checked' or validation['states']!=10
            or validation['original_project_installed'] or validation['quantum_calculations_run']
            or validation['full_historical_source_closure_revalidated'] or validation['public_archive_deposited']
            or validation['scientific_dependency']!='numpy==2.2.6' or len(validation['blocked_original_paths'])!=4
            or sha(package/'manifest.json')!=validation['package_manifest_sha256']
            or sha(folder/'exact_molecular_CPU_manifest.json')!=validation['package_manifest_sha256']
            or sha(package/'requirements.lock.txt')!=validation['scientific_requirements_sha256']):
        raise ValueError('Frozen isolated CPU example identity/scope differs')
    script=root/'examples/exact_viologen_replay/replay.py'
    if sha(script)!=manifest['files_sha256']['replay.py']:raise ValueError('Maintained and frozen example runners differ')
    definition=importlib.util.spec_from_file_location('checked_exact_molecular_CPU_replay',script)
    reader=importlib.util.module_from_spec(definition);definition.loader.exec_module(reader);current=reader.replay(package)
    if (current['states']!=10 or current['Hessians_reconstructed']!=9 or current['FD_states_reconstructed']!=6
            or current['material_admission'] or current['global_conformer_convergence_proven']
            or current['selected_precursor_id']!=frozen['selected_precursor_id']):raise ValueError('CPU replay coverage or selection differs')
    for cutoff in ('25','50','100'):
        if set(current['by_cutoff'][cutoff])!=set(frozen['by_cutoff'][cutoff]):raise ValueError('Cached observable coverage differs')
        for key,value in current['by_cutoff'][cutoff].items():
            expected=frozen['by_cutoff'][cutoff][key]
            if not math.isfinite(value) or not math.isfinite(expected) or abs(value-expected)>1e-6:raise ValueError('Cached molecular quantity differs')
    return dict(status='ten_state_exact_molecular_CPU_example_reconstructed',states=10,Hessians=9,FD_states=6,
        first_reduction_V=current['by_cutoff']['50']['E1_V_vs_Fc'],grafting_kJ_mol=current['by_cutoff']['50']['graft_precursor_2_kJ_mol'],
        original_paths_blocked=4,full_historical_source_closure_revalidated=False,public_archival_DOI=False)
