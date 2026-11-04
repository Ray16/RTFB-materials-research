"""Check the frozen fresh xTB-component method evidence without quantum imports."""
import json
import math
from pathlib import Path


def audit(directory):
    directory=Path(directory);run=json.loads((directory/'xtb_runtime_installation.json').read_text())
    inventory=json.loads((directory/'xtb_runtime_packages.json').read_text());installed=json.loads((directory/'xtb_installed_packages.json').read_text())
    numeric=run['numerical_checks'];prefix=Path(run['qualified_runtime']['prefix'])
    expected=dict(G_Ha=.017417501635698,H_Ha=.050275771916811,TS_Ha=.032858270281112)
    packages={row['name']:row for row in inventory['packages']}
    if (run['status']!='fresh_locked_xtb_runtime_and_relocated_thermo_checked' or not run['original_xtb_files_masked']
            or run['installed_packages']!=18 or len(packages)!=18 or set(packages)!=set(installed)
            or run['Python_environment_recreated'] or run['full_quantum_environment_qualified'] or numeric['CUDA_imported']
            or (packages['xtb']['version'],packages['xtb']['build'])!=('6.7.1','gfortran_ha02e54d_5')
            or any(installed[name]['sha256']!=row['sha256'] for name,row in packages.items())
            or any(not math.isfinite(numeric['upstream_kernel'][key]) or abs(numeric['upstream_kernel'][key]-value)>1e-12 for key,value in expected.items())
            or not numeric['soft_mode_retained'] or abs(numeric['soft_frequency_cm']-.391)>1e-10
            or numeric['library']!=str(prefix/'lib/libxtb.so.6.7.1') or numeric['executable']!=str(prefix/'bin/xtb')
            or run['qualified_runtime']['library_sha256']!=inventory['runtime']['library_sha256']
            or run['qualified_runtime']['binary_sha256']!=inventory['runtime']['binary_sha256']):
        raise ValueError('Frozen fresh xTB component evidence does not support the stated scope')
    return dict(status='fresh_locked_xtb_component_evidence_reconstructed',packages=18,upstream_kernel_tolerance_Ha=1e-12,
        retained_soft_frequency_cm=numeric['soft_frequency_cm'],original_paths_masked=True,full_quantum_environment_qualified=False)
