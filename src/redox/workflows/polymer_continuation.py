"""Source-preserving continuations of finite, neutral polymer calculations."""

from contextlib import contextmanager
import gc
import math


def review_memory_pilot(report, result, device_mib, reserve_mib=4096):
    """Review completed numerical work independently of the prior resource guard.

    This does not turn an unconverged energy or geometry into a qualified result.
    The input geometry is a resource control, never a physical minimum.
    """
    if (
        report.get("error")
        != "RuntimeError('Large-model memory pilot exceeds admitted safety budget')"
    ):
        raise ValueError(
            "Only a completed calculation stopped by this guard is reusable"
        )
    if (
        result["result"].get("converged_smd") is not True
        or result["result"].get("converged_gas") is not True
    ):
        raise ValueError("Pilot SCFs did not converge")
    for value in (
        result["result"].get("e_smd_eV"),
        result["result"].get("e_gas_eV"),
        result.get("gradient_max_Ha_Bohr"),
    ):
        if type(value) not in (float, int) or not math.isfinite(value):
            raise ValueError("Pilot numerical evidence is incomplete")
    peak = report.get("observed_peak_GPU_MiB")
    if (
        type(peak) is not int
        or peak <= 0
        or reserve_mib < 4096
        or peak + reserve_mib > device_mib
    ):
        raise ValueError("Insufficient measured GPU memory headroom")
    return dict(
        numerical_pilot_completed=True,
        observed_peak_GPU_MiB=peak,
        device_memory_MiB=device_mib,
        available_headroom_MiB=device_mib - peak,
        required_headroom_MiB=reserve_mib,
        maximum_memory_pilot_MiB=device_mib - reserve_mib,
        minimum_device_MiB=device_mib,
        qualification_scope=(
            "Neutral energy/gradient execution only; no minimum or Hessian claim"
        ),
    )


def separate_stages(original, *args, **kwargs):
    """Release the optimizer's SCF/scanner before building the larger SP object.

    The old monolithic call kept its geometry mean field alive during the larger
    active SP. The two calls here use identical equations, solver settings and
    final coordinates, and retain ordinary SP adoption/provenance.
    """
    if not kwargs.get("do_opt", True):
        return original(*args, **kwargs)
    if kwargs.get("do_freq", True) or not kwargs.get("opt_out"):
        raise ValueError(
            "This adapter requires explicit opt_out and separate thermal evaluation"
        )
    geometry_options = dict(kwargs, do_gas=False, do_smd=False, do_freq=False)
    geometry = original(*args, **geometry_options)
    if geometry.get("optimized") is not True:
        raise ValueError("Optimizer did not qualify its geometry")
    gc.collect()
    from redox.qm.dft import _free_gpu_pool

    _free_gpu_pool()
    point_options = dict(kwargs, do_opt=False, do_freq=False, opt_out=None)
    point = original(kwargs["opt_out"], *args[1:], **point_options)
    merged = dict(geometry, **point)
    for field in (
        "xyz",
        "optimized",
        "geometry_optimization_converged",
        "optimization_SCF_converged",
        "opt_gradient_grid_response",
        "opt_out",
    ):
        if field in geometry:
            merged[field] = geometry[field]
    merged["execution_stages_separated"] = True
    return merged


@contextmanager
def execution_settings(optimizer=None):
    """Run-local hooks; never edit a producer used by another live campaign."""
    import redox.qm.dft as dft

    energy, builder, driver = dft.dft_smd, dft._build_smd_mf, dft._geomopt_fn

    def evaluate(*args, **kwargs):
        return separate_stages(energy, *args, **kwargs)

    def build(*args, **kwargs):
        mf = builder(*args, **kwargs)
        if optimizer and (args[1] if len(args) > 1 else kwargs.get("xc")) == "r2scan":
            mf.conv_tol_grad = 1e-8
        return mf

    def geometry(backend):
        kernel = driver(backend)

        def run(method, **kwargs):
            if optimizer:
                kwargs.update(coordsys="cart", trust=0.01, tmin=1e-6)
            return kernel(method, **kwargs)

        return run

    try:
        dft.dft_smd, dft._build_smd_mf, dft._geomopt_fn = evaluate, build, geometry
        yield
    finally:
        dft.dft_smd, dft._build_smd_mf, dft._geomopt_fn = energy, builder, driver
