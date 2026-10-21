"""Guard reuse requires numerical success; optimizer and SP lifetimes are separate."""

from copy import deepcopy
import pytest

from redox.workflows.polymer_continuation import review_memory_pilot, separate_stages


def evidence():
    report = dict(
        error="RuntimeError('Large-model memory pilot exceeds admitted safety budget')",
        observed_peak_GPU_MiB=89062,
    )
    result = dict(
        result=dict(
            converged_smd=True, converged_gas=True, e_smd_eV=-100.0, e_gas_eV=-99.0
        ),
        gradient_max_Ha_Bohr=0.01,
    )
    return report, result


def test_guard_only_stop_can_be_reviewed_with_actual_device_headroom():
    report, result = evidence()
    reviewed = review_memory_pilot(report, result, 95830)
    assert reviewed["available_headroom_MiB"] == 6768
    assert reviewed["maximum_memory_pilot_MiB"] == 91734
    assert "no minimum" in reviewed["qualification_scope"]


@pytest.mark.parametrize("problem", ["SCF", "gradient", "other_failure", "memory"])
def test_review_never_promotes_incomplete_work_or_insufficient_memory(problem):
    report, result = evidence()
    device = 95830
    if problem == "SCF":
        result["result"]["converged_smd"] = False
    elif problem == "gradient":
        result["gradient_max_Ha_Bohr"] = float("nan")
    elif problem == "other_failure":
        report["error"] = "SCF failed"
    else:
        device = 90000
    with pytest.raises(ValueError):
        review_memory_pilot(report, result, device)


def test_separate_calls_retain_exact_opt_geometry_and_active_SP_provenance(
    monkeypatch, tmp_path
):
    import redox.qm.dft as dft

    released = []
    monkeypatch.setattr(dft, "_free_gpu_pool", lambda: released.append(True))
    calls = []

    def original(path, charge, mult, **kwargs):
        calls.append((path, charge, mult, deepcopy(kwargs)))
        if kwargs["do_opt"]:
            return dict(
                optimized=True,
                xyz=str(path),
                opt_out=str(kwargs["opt_out"]),
                geometry_optimization_converged=True,
                optimization_SCF_converged=True,
            )
        assert released
        return dict(
            optimized=False,
            e_smd_eV=-10.0,
            e_gas_eV=-9.0,
            energy_geom_sha1="final_geometry",
        )

    initial, optimized = tmp_path / "input.xyz", tmp_path / "opt.xyz"
    result = separate_stages(
        original,
        initial,
        0,
        1,
        do_opt=True,
        do_freq=False,
        opt_out=optimized,
        sp_basis="def2-tzvpd",
    )
    assert (
        len(calls) == 2
        and calls[0][3]["do_smd"] is False
        and calls[0][3]["do_gas"] is False
    )
    assert calls[1][0] == optimized and calls[1][3]["sp_basis"] == "def2-tzvpd"
    assert (
        result["optimized"] is True and result["energy_geom_sha1"] == "final_geometry"
    )
    assert result["e_smd_eV"] == -10.0 and result["geometry_optimization_converged"]


def test_nonoptimized_resource_pilot_does_not_change_execution():
    calls = []

    def original(*args, **kwargs):
        calls.append((args, kwargs))
        return {"marker": "same call"}

    result = separate_stages(original, "input.xyz", 0, 1, do_opt=False, do_freq=False)
    assert result == {"marker": "same call"} and len(calls) == 1
