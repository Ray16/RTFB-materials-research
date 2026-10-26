import importlib.util
from pathlib import Path

import pytest
from pyscf import df, gto

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('fixed_counterpoise_aux_test',ROOT/'scripts/validation/run_cpu_counterpoise_reference.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_actual_df_build_and_range_separation_use_complete_auxiliary_space():
    own=gto.M(atom='H 0 0 0; H 0 0 .74',basis='def2-svp',verbose=0)
    full=gto.M(atom='H 0 0 0; H 0 0 .74; H 0 0 3; H 0 0 3.74',basis='def2-svp',verbose=0)
    outsider=gto.M(atom='He 0 0 0',basis='def2-svp',verbose=0)
    original=df.addons.make_auxmol
    full_n=original(full,None).nao_nr();own_n=original(own,None).nao_nr()
    assert full_n>own_n
    with module.fixed_complex_auxiliary_factory(df.addons,full,[full,own]):
        fit=df.DF(own);fit.build()
        assert fit.auxmol.nao_nr()==full_n
        with fit.range_coulomb(.3) as long_range:
            long_range.build()
            assert long_range.auxmol.nao_nr()==full_n
            assert long_range.auxmol.omega==pytest.approx(.3)
        assert df.addons.make_auxmol(outsider,None).nao_nr()==original(outsider,None).nao_nr()
        with pytest.raises(ValueError):df.addons.make_auxmol(own,'def2-universal-jkfit')
    assert df.addons.make_auxmol is original
    restored=df.DF(own);restored.build();assert restored.auxmol.nao_nr()==own_n


def test_factory_restores_on_failure_and_rejects_wrong_complex():
    one=gto.M(atom='He 0 0 0',basis='def2-svp',verbose=0)
    other=gto.M(atom='Ne 0 0 0',basis='def2-svp',verbose=0)
    original=df.addons.make_auxmol
    with pytest.raises(ValueError):
        with module.fixed_complex_auxiliary_factory(df.addons,one,[other]):pass
    with pytest.raises(RuntimeError):
        with module.fixed_complex_auxiliary_factory(df.addons,one,[one]):raise RuntimeError('test')
    assert df.addons.make_auxmol is original
