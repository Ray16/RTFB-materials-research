import copy
import json

import pytest
from rdkit import Chem

from redox.build.contracts import validate_group
from redox.build.structures import decorate
from redox.core import common, protocol
from redox.screening.pareto import _objectives
from redox.validation.readiness import inspect_states

SCAFFOLD = "Cc1ccc(C[*:1])cc1"
MOTIF = dict(id="new_quinone", name="test quinone", family="quinone (n-type)",
             frag="[*:1]OC1=CC(=O)C=CC1=O",
             states=[("neu", 0, 1, 0), ("red1", -1, 2, -1), ("red2", -2, 1, -2)])


def test_new_motif_atom_balance_and_states():
    mol = validate_group(SCAFFOLD, MOTIF)
    expected = Chem.MolFromSmiles("Cc1ccc(COC2=CC(=O)C=CC2=O)cc1")
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(expected)


@pytest.mark.parametrize("frag", ["OC", "[*:2]OC", "[*:1]O.C", "[*:1]OC[*:2]"])
def test_bad_attachment_is_rejected(frag):
    with pytest.raises(ValueError):
        decorate(SCAFFOLD, frag)


@pytest.mark.parametrize("change", [
    {"id": "../bad"}, {"states": [("neu", 0, 2, 0), ("red", -1, 1, -1)]},
    {"states": [("neu", 0, 1, 0), ("red", -1, 2, 1)]},
    {"states": [("neu", 0, 1, 0), ("red", -2, 1, -2)]},
    {"frag": "[*:1]OC(F)Cl"},
])
def test_invalid_motif_fails_before_build(change):
    group = copy.deepcopy(MOTIF); group.update(change)
    with pytest.raises(ValueError):
        validate_group(SCAFFOLD, group)


def test_never_computed_state_is_incomplete(tmp_path):
    rows = inspect_states([dict(id="new", state="neu", charge=0)], root=tmp_path)
    assert rows[0]["status"] == "INCOMPLETE"
    assert "missing_geometry" in rows[0]["issues"]


def test_readiness_uses_the_actual_active_energy_for_bound_ensemble_policy(tmp_path,monkeypatch):
    from redox.validation import readiness
    from redox.qm import dft
    sd=tmp_path/'m/neu';sd.mkdir(parents=True)
    (sd/'opt.xyz').write_text('2\n\nH 0 0 0\nH 0 0 .74\n')
    (sd/'result.json').write_text(json.dumps(dict(charge=0,mult=1)))
    monkeypatch.setattr(dft,'_opt_protocol_mismatch',lambda *args:'')
    monkeypatch.setattr(readiness,'load_record',lambda *args:dict(e_smd_eV=-11.))
    observed=[]
    def selected(record,directory):
        observed.append(record['e_smd_eV'])
        return dict(record,g_thermal_eV=1. if record['e_smd_eV']==-11. else None,n_imag=0,
                    thermal_definition='native_solution_dft_qrrho_spin',
                    thermal_qc='ok' if record['e_smd_eV']==-11. else 'stale_or_invalid_spin_ensemble')
    monkeypatch.setattr(readiness,'apply_thermal_policy',selected)
    row,=readiness.inspect_states([dict(id='m',state='neu',charge=0)],root=tmp_path)
    assert row['status']=='ready' and observed==[-11.]


def test_partial_lambda_cannot_be_a_complete_objective():
    sc = common.load_config("scorecard_config")
    row = dict(E_V=-1, sigma_E_V=.2, specific_capacity_mAh_g=100,
               disp_applicable=False, lambda_i_eV=.1, lambda_qc="partial")
    obj, missing = _objectives(row, "anolyte", sc)
    assert "kinetics" in missing and "kinetics" not in obj


def test_nonfinite_energies_never_become_values():
    assert common.free_energy(dict(e_smd_eV=float("nan"), g_thermal_eV=1))[0] is None
    assert common.free_energy(dict(e_smd_eV=1, g_thermal_eV=float("inf")))[0] is None
    assert common.to_float("NaN") is None
    assert common.free_energy(dict(e_smd_eV=None, g_thermal_eV=1, n_imag=2))[1] == "imag_modes"


@pytest.mark.parametrize("field,value", [("charge", 1), ("geom_sha1", "wrong"),
                                        ("protocol", {}), ("e_smd_eV", float("nan")),
                                        ("converged_smd", "False")])
def test_record_contents_must_match_address(tmp_path, field, value):
    (tmp_path / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 0.7\n")
    path = protocol.record_path(tmp_path, 0, 1)
    path.parent.mkdir()
    rec = dict(protocol=protocol.ACTIVE_SP, protocol_hash=path.stem,
               geom_sha1=protocol.geom_sha1(tmp_path / "opt.xyz"), charge=0, mult=1,
               e_smd_eV=-1, e_gas_eV=-1, converged_smd=True, converged_gas=True)
    path.write_text(json.dumps(rec))
    assert protocol.load_record(tmp_path, 0, 1)
    rec[field] = value; path.write_text(json.dumps(rec))
    assert protocol.load_record(tmp_path, 0, 1) is None


def test_builder_real_embedding_resume_and_collision(tmp_path, monkeypatch):
    import sys
    from redox.build import candidates
    config = tmp_path / "motifs.py"
    config.write_text(f"SCAFFOLD = {SCAFFOLD!r}\nGROUPS = {[MOTIF]!r}\n")
    library = tmp_path / "library"
    monkeypatch.setattr(sys, "argv", ["build", "--config", str(config), "--library", str(library)])
    candidates.main()
    before = {p: p.read_bytes() for p in library.rglob("*") if p.is_file()}
    assert len(list((library / MOTIF["id"] / "conformers").glob("*.xyz"))) > 0
    candidates.main()
    assert before == {p: p.read_bytes() for p in library.rglob("*") if p.is_file()}
    changed = dict(MOTIF, frag="[*:1]OC1=C(F)C(=O)C=CC1=O")
    config.write_text(f"SCAFFOLD = {SCAFFOLD!r}\nGROUPS = {[changed]!r}\n")
    with pytest.raises(ValueError, match="different chemistry"):
        candidates.main()
    assert before == {p: p.read_bytes() for p in library.rglob("*") if p.is_file()}


def test_ambiguous_primary_source_is_not_tier_a_spacing():
    rows = common.load_config("benchmark_mecn").BENCHMARK
    entry = next(r for r in rows if r["name"] == "2-fluoro-1,4-benzoquinone")
    events = [e for e in entry["events"] if "0c03632" in e["doi"]]
    assert next(e for e in events if e["event"] == "neu->red1")["tier"] == "A"
    assert next(e for e in events if e["event"] == "red1->red2")["tier"] == "B"
    assert next(e for e in events if e["kind"] == "spacing")["tier"] == "B"


def test_multielectron_summary_does_not_add_oxidation_to_reduction():
    from redox.properties.potentials import _multielectron_summary
    events = [dict(id="radical", name="test", family="nitroxide", q_ox=1,
                   event="ox->rad", E_vs_Fc_V=.5),
              dict(id="radical", name="test", family="nitroxide", q_ox=0,
                   event="rad->red", E_vs_Fc_V=-1)]
    assert _multielectron_summary(events, (-3, 2), {"radical": 0}) == []
    events = [dict(id="q", name="quinone", family="quinone", q_ox=0,
                   event="neu->red1", E_vs_Fc_V=-1),
              dict(id="q", name="quinone", family="quinone", q_ox=-1,
                   event="red1->red2", E_vs_Fc_V=-2)]
    summary = _multielectron_summary(events, (-3, 2), {"q": 0})
    assert len(summary) == 1 and summary[0]["usable_electrons"] == 2


def test_failed_sp_cache_is_preserved_and_recomputed(tmp_path, monkeypatch):
    from redox.qm import sp, dft
    (tmp_path / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 0.7\n")
    (tmp_path / "result.json").write_text(json.dumps(dict(charge=0, mult=1)))
    path = protocol.record_path(tmp_path, 0, 1)
    path.parent.mkdir()
    old = dict(protocol=protocol.ACTIVE_SP, protocol_hash=path.stem,
               geom_sha1=protocol.geom_sha1(tmp_path / "opt.xyz"), charge=0, mult=1,
               converged_smd=False, converged_gas=True, e_smd_eV=-1, e_gas_eV=-1)
    path.write_text(json.dumps(old))
    monkeypatch.setattr(dft, "dft_smd", lambda *args, **kwargs: dict(
        converged_smd=True, converged_gas=True, e_smd_eV=-2, e_gas_eV=-1.5))
    monkeypatch.setattr(sp, "software_versions", lambda: {})
    status, result = sp.ensure("m", "neu", tmp_path, 0, 1, backend="cpu")
    assert status == "computed" and result["e_smd_eV"] == -2
    archived = list((path.parent / "_superseded").glob("*.json"))
    assert len(archived) == 1 and json.loads(archived[0].read_text()) == old
    assert protocol.load_record(tmp_path, 0, 1)["e_smd_eV"] == -2


def test_dft_does_not_report_incomplete_cache_as_success(tmp_path, monkeypatch):
    from redox.qm import dft
    sd = tmp_path / "calcs/dft/m/neu"
    sd.mkdir(parents=True)
    (sd / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 0.7\n")
    (sd / "result.json").write_text(json.dumps(dict(charge=0, mult=1,
                                                    converged_smd=False, g_thermal_eV=1)))
    monkeypatch.setattr(dft, "ROOT", tmp_path)
    monkeypatch.setattr(dft, "_opt_protocol_mismatch", lambda *args: "")
    assert dft.run_batch([dict(id="m", state="neu", charge="0")], force=False) == 1


@pytest.mark.parametrize("change", [{"e_gas_eV": float("nan")}, {"e_gas_eV": float("inf")},
                                   {"protocol": dict(protocol.ACTIVE_SP, density_fit=False)},
                                   {"protocol": dict(protocol.ACTIVE_SP, conv_tol=1e-6)}])
def test_cross_point_must_have_finite_energy_and_full_provenance(change):
    from redox.properties.reorg import _cache_ok, CACHE_SCHEMA
    record = dict(e_gas_eV=-1, sp_xc="wb97m-v", sp_basis="def2-tzvpd", sp_nlc="vv10",
                  sp_disp=None, cache_schema=CACHE_SCHEMA, gas_homo_eV=-1,
                  anion_unbound=False, charge=-1, mult=2, geom_sha1="g",
                  converged_gas=True, protocol=protocol.ACTIVE_SP)
    assert _cache_ok(record, -1, 2, "g")
    record.update(change)
    assert not _cache_ok(record, -1, 2, "g")


def test_conformer_constraints_follow_stored_atom_order(tmp_path):
    import importlib.util
    from rdkit.Chem import AllChem
    path = common.ROOT / "scripts/analysis/reorg_conformer_matched.py"
    spec = importlib.util.spec_from_file_location("cf_mapping_test", path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    # Embedding order intentionally differs from canonical SMILES.
    mol = Chem.AddHs(Chem.MolFromSmiles("OCCc1ccc(C)cc1"))
    assert AllChem.EmbedMolecule(mol, randomSeed=17) == 0
    xyz = tmp_path / "stored.xyz"
    xyz.write_text(Chem.MolToXYZBlock(mol))
    output = tmp_path / "freeze.txt"
    n = mod.rotatable_dihedral_constraints(Chem.MolToSmiles(Chem.RemoveHs(mol)), output, xyz)
    assert n > 0
    for line in output.read_text().splitlines()[1:]:
        quad = [int(i) - 1 for i in line.split()[1:]]
        assert all(mol.GetAtomWithIdx(i).GetAtomicNum() > 1 for i in quad)
        assert all(mol.GetBondBetweenAtoms(i, j) is not None for i, j in zip(quad, quad[1:]))


def test_geometry_skeleton_mapping_ignores_radical_electron_labels(tmp_path):
    from rdkit.Chem import AllChem
    molecule = Chem.AddHs(Chem.MolFromSmiles("CC[O]"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=19) == 0
    xyz = tmp_path / "radical.xyz"; xyz.write_text(Chem.MolToXYZBlock(molecule))
    mapping = common.smiles_to_xyz_map("CC[O]", xyz)
    assert mapping is not None and len(mapping) == 3
    skeleton = common._skeleton(Chem.MolFromSmiles("CC[O]"))
    assert all(atom.GetNumRadicalElectrons() == 0 for atom in skeleton.GetAtoms())


def test_thermal_cache_follows_geometry_and_reference_changes(tmp_path):
    import importlib.util
    path = common.ROOT / "scripts/analysis/thermal_sensitivity.py"
    spec = importlib.util.spec_from_file_location("thermal_contract_test", path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    xyz = tmp_path / "opt.xyz"; xyz.write_text("2\n\nH 0 0 0\nH 0 0 0.7\n")
    cache, key = mod.thermal_cache_path(tmp_path, 0, 1, 50)
    assert cache != mod.thermal_cache_path(tmp_path, 0, 1, 25)[0]
    xyz.write_text("2\n\nH 0 0 0\nH 0 0 0.8\n")
    assert cache != mod.thermal_cache_path(tmp_path, 0, 1, 50)[0]
    # Equal shifts of molecular and Fc redox free energies cancel on the Fc scale.
    base = mod.referenced_potential(-95, -100, -45, -50)
    shifted = mod.referenced_potential(-94.8, -100, -44.8, -50)
    assert base == pytest.approx(shifted)
    assert mod.referenced_potential(None, -100, -45, -50) is None


def test_ready_energy_does_not_hide_a_changed_declared_structure(tmp_path, monkeypatch):
    from rdkit.Chem import AllChem
    from redox.qm import dft
    mol = Chem.AddHs(Chem.MolFromSmiles("COC"))
    assert AllChem.EmbedMolecule(mol, randomSeed=29) == 0
    sd = tmp_path / "m/neu"; sd.mkdir(parents=True)
    (sd / "opt.xyz").write_text(Chem.MolToXYZBlock(mol))
    (sd / "result.json").write_text(json.dumps(dict(charge=0, mult=1, g_thermal_eV=1)))
    p = protocol.record_path(sd, 0, 1); p.parent.mkdir()
    p.write_text(json.dumps(dict(protocol=protocol.ACTIVE_SP, protocol_hash=p.stem,
                                geom_sha1=protocol.geom_sha1(sd / "opt.xyz"), charge=0, mult=1,
                                converged_smd=True, converged_gas=True, e_smd_eV=-2, e_gas_eV=-1)))
    monkeypatch.setattr(dft, "_opt_protocol_mismatch", lambda *args: "")
    result = inspect_states([dict(id="m", state="neu", charge=0, smiles="CCO")], root=tmp_path)[0]
    assert result["status"] == "INVALID" and "declared_connectivity_mismatch" in result["issues"]


def test_spin_gap_uses_active_level_and_does_not_validate_imaginary_thermal(tmp_path, monkeypatch):
    from redox.screening import scorecard
    from redox.validation import spins
    from redox.qm import dft
    helper = tmp_path / "spincheck/m/red2_m3"; helper.mkdir(parents=True)
    (helper / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 .8\n")
    (helper / "result.json").write_text(json.dumps(dict(charge=-2, mult=3, e_smd_eV=-50,
                                                       g_thermal_eV=7, n_imag=2)))
    address = protocol.record_path(helper, -2, 3); address.parent.mkdir()
    address.write_text(json.dumps(dict(protocol=protocol.ACTIVE_SP, protocol_hash=address.stem,
                                       geom_sha1=protocol.geom_sha1(helper / "opt.xyz"), charge=-2, mult=3,
                                       converged_gas=True, converged_smd=True, e_gas_eV=-98, e_smd_eV=-99)))
    monkeypatch.setattr(scorecard, "SPINCHECK", tmp_path / "spincheck")
    monkeypatch.setattr(scorecard, "UMA", tmp_path / "uma")
    monkeypatch.setattr(common, "read_result", lambda *a, **k: dict(charge=-2, mult=1, e_smd_eV=-100,
                                                                  g_thermal_eV=20, n_imag=2))
    monkeypatch.setattr(spins, 'snapshot', lambda *args: {'geometry': 'current'})
    monkeypatch.setattr(spins, 'read_comparison', lambda *args: dict(status='spin_gap_verified',
        chosen_identity={'geometry': 'current'}, alternative_mult=3, electronic_gap_eV=1))
    monkeypatch.setattr(dft, '_opt_protocol_mismatch', lambda *args: '')
    gap, source = scorecard._min_spin_gap_eV("m")
    assert gap == 1 and source == "dft_electronic(thermal_unvalidated):red2"
    monkeypatch.setattr(common, "read_result", lambda *a, **k: dict(charge=-2,mult=3,e_smd_eV=-99,
                                                                  g_thermal_eV=None,n_imag=0))
    assert scorecard._min_spin_gap_eV('m')==(None,'')  # selected state is not its own alternative
    address.unlink()
    assert scorecard._min_spin_gap_eV("m") == (None, "")


def test_occupation_seed_can_find_lower_solution_than_standard_guesses(monkeypatch):
    from redox.qm import dft
    monkeypatch.setattr(dft, "_free_gpu_pool", lambda: None)
    class FakeMF:
        pass
    def kernel(mf, dm0=None):
        return (-3.0, True) if dm0 == 'ground_state_density' else (-1.0, True)
    monkeypatch.setattr(dft, '_kernel_robust', kernel)
    e, converged, _, evidence = dft._lowest_scf(FakeMF, ('minao','atom'),
                                               initial_densities={'rm0':'ground_state_density'})
    assert e == -3 and converged
    assert [r['guess'] for r in evidence] == ['minao','atom','rm0']


def test_newton_fallback_energy_and_orbitals_stay_on_same_solution():
    from redox.qm.dft import _kernel_robust
    class Newton:
        mo_coeff='new_density';mo_occ='new_occ';mo_energy='new_levels';e_tot=-3;converged=True
        def kernel(self,dm0=None): return self.e_tot
    class MF:
        converged=False
        def kernel(self,dm0=None): return -1
        def make_rdm1(self): return 'old_density'
        def newton(self): return Newton()
    mf=MF()
    e,converged=_kernel_robust(mf)
    assert e == mf.e_tot == -3 and converged and mf.converged
    assert mf.mo_coeff=='new_density'


def test_unvalidated_thermal_terms_cannot_define_complete_pareto_objectives():
    from redox.screening.pareto import _objectives
    sc=common.load_config('scorecard_config')
    candidate=dict(E_V=-1, sigma_E_V=.2, specific_capacity_mAh_g=100, dG_disp_kJmol=60,
                   disp_applicable=True, lambda_i_eV=.5, lambda_qc='ok', thermal_qc='imag_modes')
    objectives,missing=_objectives(candidate,'anolyte',sc)
    assert {'voltage','stability'}.issubset(missing)
    assert 'voltage' not in objectives and 'stability' not in objectives
    assert 'capacity' in objectives and 'kinetics' in objectives


def test_xtb_hessian_writer_matches_fixed_width_column_reader(tmp_path):
    import numpy as np
    from redox.qm.vibrations import write_xtb_hessian
    matrix=np.arange(36,dtype=float).reshape(6,6)/10
    out=tmp_path/'hessian';write_xtb_hessian(out,matrix)
    lines=out.read_text().splitlines()[1:]
    assert len(lines)==12
    columns=[]
    for index in range(0,12,2):
        values=[]
        for line in lines[index:index+2]:
            assert line.startswith('     ')
            values.extend(float(line[start:start+15]) for start in range(5,len(line),15))
        columns.append(values)
    np.testing.assert_allclose(np.array(columns).T,matrix)


def test_thermal_parser_uses_full_precision_rrho_term():
    from redox.qm.vibrations import parse_xtb_thermal_eV
    text='298.15 .02 .328722 .067 .261548\n :: G(RRHO) contrib. 0.261547738424 Eh ::'
    assert parse_xtb_thermal_eV(text)==pytest.approx(.261547738424*common.HARTREE_EV,abs=1e-12)


def test_uma_gas_diagnostic_uses_same_model_and_same_phase_fc(monkeypatch):
    from redox.properties import potentials
    energies={('m','ox'):-20,('m','red'):-25,('ferrocene','ox'):-10,('ferrocene','neu'):-16}
    charges={('m','ox'):0,('m','red'):-1,('ferrocene','ox'):1,('ferrocene','neu'):0}
    def read(gid,st,root=None):
        return dict(energy_eV=energies[gid,st],charge=charges[gid,st],converged=True,model='same_model')
    monkeypatch.setattr(potentials,'read_result',read)
    delta,value,qc=potentials.uma_gas_diagnostic('m','ox','red',0,-1)
    assert delta==-5 and value==-1 and 'not_MeCN_G' in qc
    assert potentials.uma_gas_diagnostic('m','ox','red',1,0)[0] is None


def test_neutral_solvation_proxy_never_substitutes_a_charged_species():
    from redox.properties.capacity_and_proxies import _dG_solv_neutral
    assert _dG_solv_neutral([('ox',2,{'dG_solv_eV':-5}),('red',1,{'dG_solv_eV':-2})]) is None
    assert _dG_solv_neutral([('neu',0,{'dG_solv_eV':-.5}),('red',-1,{'dG_solv_eV':-2})])==-.5


def test_unknown_small_spin_gap_excludes_a_complete_primary_front():
    from redox.screening.pareto import _objectives
    row=dict(E_V=-1,sigma_E_V=.2,specific_capacity_mAh_g=100,lambda_i_eV=.4,lambda_qc='ok',
             disp_applicable=False,thermal_qc='ok',spin_confidence='DFT_check_required')
    _,missing=_objectives(row,'anolyte',common.load_config('scorecard_config'))
    assert 'spin_state' in missing


def test_verified_dft_spin_decision_overrides_uma_and_rejects_changed_species(tmp_path,monkeypatch):
    import hashlib
    from redox.qm import dft
    row=dict(id='m',state='red2',charge='-2',smiles='CC')
    p=tmp_path/'config/spin_overrides.json';p.parent.mkdir()
    p.write_text(json.dumps({'m':{'red2':dict(mult=1,charge=-2,
                                             smiles_sha256=hashlib.sha256(b'CC').hexdigest(),
                                             protocol=protocol.ACTIVE_SP)}}))
    uma=tmp_path/'calcs/uma/m/red2/result.json';uma.parent.mkdir(parents=True)
    uma.write_text(json.dumps(dict(chosen_mult=3)))
    monkeypatch.setattr(dft,'ROOT',tmp_path);monkeypatch.setattr(dft,'read_manifest',lambda:[row])
    assert dft._uma_mult('m','red2',3)==1
    row['smiles']='CO'
    with pytest.raises(RuntimeError,match='revalidation'):dft._uma_mult('m','red2',3)


def test_unstable_rrho_produces_only_disproportionation_estimates(monkeypatch):
    from redox.properties import stability
    energies={'ox':0.,'mid':-2.,'red':-3.}
    monkeypatch.setattr(stability,'_states_by_charge',lambda gid:[('ox',0),('mid',-1),('red',-2)])
    monkeypatch.setattr(stability,'read_result',lambda gid,st,**kw:dict(e_smd_eV=energies[st],
                        g_thermal_eV=0,n_imag=2 if st=='red' else 0))
    row=stability.disproportionation('m')[0]
    assert row['status']=='THERMAL_UNVALIDATED' and 'dG_disp_eV' not in row
    assert row['dG_disp_eV_estimate']==1
