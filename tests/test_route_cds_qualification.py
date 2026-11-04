import json
import numpy as np
import pytest
from redox.core.hashing import sha256_file
from redox.workflows.route_cds_native import check_components


@pytest.mark.parametrize('bad',[None,'fine_nan','wrong_shape','gradient_nan'])
def test_all_component_arrays_and_gradient_checks_must_be_finite(tmp_path,bad):
    paths=[]
    for sign in ('minus','plus'):
        path=tmp_path/(sign+'.json');arrays=path.with_suffix('.npz');base=np.zeros((117,117));base[0,1]=1e-5
        fine=np.zeros((117,117))
        if sign=='minus' and bad=='fine_nan':fine[0,0]=np.nan
        if sign=='minus' and bad=='wrong_shape':fine=np.zeros((1,1))
        np.savez(arrays,saved_full=base,residual_without_default_CDS=np.zeros((117,117)),
            **{'CDS_0.0001':base,'CDS_1e-05':np.zeros((117,117)),'CDS_1e-06':fine,'CDS_1e-07':np.zeros((117,117))})
        error=float('nan') if sign=='minus' and bad=='gradient_nan' else 0.
        row=dict(status='legacy_CDS_Hessian_step_component_study_not_full_QC',key='phenothiazine_neu_mode_'+sign,
            native_G_qualified=False,CUDA_initialized=False,CDS_energy_and_gradient_functions_changed=False,
            electrostatic_model_changed=False,source_sha256={str(arrays):sha256_file(arrays)},
            energy_gradient_checks=[dict(tested_components=117,maximum_energy_gradient_error_Ha_Bohr=error)]*2,
            geometry_sha256='fixture',fine_pair_maximum_change_Ha_Bohr2=0.,
            residual_without_default_CDS_maximum_asymmetry_Ha_Bohr2=0.)
        path.write_text(json.dumps(row));paths.append(path)
    if bad:
        with pytest.raises(ValueError):check_components(paths)
    else:
        result=check_components(paths)
        assert result['step_Bohr']==1e-6 and not result['native_G_qualified']
