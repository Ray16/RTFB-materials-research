"""Check the versioned nine-gradient-panel corrected molecular CPU cache."""
import json
import types


def audit(root,evidence):
    package=root/'examples/exact_viologen_corrected_replay/data';source=package/'replay.py'
    runner=types.ModuleType('checked_corrected_example');runner.__file__=str(source)
    exec(compile(source.read_bytes(),str(source),'exec'),runner.__dict__)
    result=runner.replay(package)
    if (result['states'],result['Hessians_reconstructed'],result['FD_states_reconstructed'])!=(10,9,9):
        raise ValueError('Corrected molecular cache coverage differs')
    cycles=json.loads((evidence/'primary_native_cohort_cycles.json').read_text())
    for row in cycles['cycles']:
        for cutoff,value in row['corrected_by_cutoff'].items():
            if abs(value-result['by_cutoff'][cutoff][row['id']])>1e-7:
                raise ValueError('Corrected CPU cache and frozen primary cohort differ')
    return dict(states=10,Hessians=9,saved_gradient_panels=9,cached_numerical_replay_passed=True,
        original_quantum_execution_repeated=False,wavefunction_spin_reconstructed=False,
        maximum_frequency_error_cm=max(r['frequency_error_cm'] for r in result['records']),
        maximum_nuclear_error_eV=max(max(r['nuclear_errors_eV'].values()) for r in result['records']))
