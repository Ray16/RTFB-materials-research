#!/usr/bin/env python
"""Replay-gated bounded PCM contractions for the checked spacing cohort."""

import argparse
from contextlib import contextmanager
import json
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
SOURCE = ROOT / "scripts/pipeline/resume_polymer_thermal.py"
MODULE = ROOT / "src/redox/qm/pcm_bounded.py"
CAMPAIGN = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_bounded_pcm"


def configure():
    base = api.load_module("bounded_polymer_base", SOURCE)
    base.PARENT = (
        ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_memory_resume"
    )
    base.CAMPAIGN, base.SCRIPT = CAMPAIGN, SCRIPT

    @contextmanager
    def solvent_memory(metrics):
        import cupy as cp
        from pyscf.solvent.hessian import pcm as cpu
        from gpu4pyscf.solvent.hessian import pcm as gpu
        from redox.qm.pcm_bounded import (
            DiagonalOperator,
            CurvatureBlocks,
            OffDiagonalTensor,
        )
        from redox.qm.hessian_memory import pcm_certificate
        from redox.qm.native_thermal import NATIVE_THERMAL

        _, certificate = pcm_certificate(NATIVE_THERMAL)
        for field, path in (
            ("cpu_pcm_library_sha256", cpu.__file__),
            ("gpu_pcm_library_sha256", gpu.__file__),
        ):
            if api.digest(path) != certificate[field]:
                raise RuntimeError("PCM primitive source changed")
        previous = (gpu.get_d2Sii, gpu.get_d2D_d2S, gpu.get_v_dot_d2A_dot_q)
        counters = dict(
            adapter_sha256=api.digest(MODULE),
            uploaded_blocks=0,
            maximum_uploaded_block_bytes=0,
            primitive_calls=0,
            peak_cached_host_bytes=0,
            calls=[],
        )
        metrics.update(
            area_calls=0,
            maximum_area_block_bytes=0,
            diagonal_contractions=0,
            host_curvature=counters,
        )

        def diagonal(surface, first, second, stream=None):
            return DiagonalOperator(surface, first, second, cp, metrics)

        def curvature(surface, with_S=True, with_D=False, stream=None):
            host = {
                name: cp.asnumpy(surface[name])
                for name in ("charge_exp", "grid_coords", "norm_vec")
            }
            n = len(host["charge_exp"])
            # Bound concurrency using cache capacity plus fixed host headroom.
            projected = min(18 * n * n * 8, 128 * 1024**3) + 128 * 1024**3
            counters["calls"].append(
                dict(
                    n_surface_points=n, projected_conservative_host_peak_bytes=projected
                )
            )
            blocks = CurvatureBlocks(host, cpu.get_d2D_d2S, cp.asarray, counters)
            return (
                OffDiagonalTensor(blocks, 0) if with_D else None
            ), OffDiagonalTensor(blocks, 1)

        def area(tensor, left, right):
            metrics["area_calls"] += 1
            metrics["maximum_area_block_bytes"] = max(
                metrics["maximum_area_block_bytes"], tensor[0].nbytes
            )
            return base.area_contraction(tensor, left, right, cp)

        gpu.get_d2Sii, gpu.get_d2D_d2S, gpu.get_v_dot_d2A_dot_q = (
            diagonal,
            curvature,
            area,
        )
        try:
            yield
        finally:
            gpu.get_d2Sii, gpu.get_d2D_d2S, gpu.get_v_dot_d2A_dot_q = previous

    base.solvent_memory = solvent_memory
    checked = base.checked_replays

    def checked_replays():
        evidence = checked()
        for path in evidence:
            metrics = json.loads(Path(path).read_text())["metrics"]
            if not metrics.get("diagonal_contractions") or not metrics[
                "host_curvature"
            ].get("primitive_calls"):
                raise ValueError(
                    "Bounded PCM derivatives were not exercised by both full replays"
                )
            if metrics["host_curvature"]["adapter_sha256"] != api.digest(MODULE):
                raise ValueError("Bounded PCM source differs from its replay")
        return evidence

    base.checked_replays = checked_replays
    original_runner = base.runner

    def runner():
        import redox.qm.projected_thermal as projected

        original = projected.ensure_projected_thermal

        def seal(path, geometry):
            record = json.loads(Path(path).read_text())
            # This is a separately replayed implementation, not the earlier
            # full-host-array adapter. Preserve truthful execution provenance.
            previous = record.pop("execution_pcm_curvature")
            record["execution_bounded_pcm"] = dict(
                adapter_sha256=api.digest(MODULE),
                driver_sha256=api.digest(SCRIPT),
                replay_evidence_sha256=checked_replays(),
                source_geometry_sha256=api.digest(geometry),
                target_evidence_sha256=previous["target_evidence_sha256"],
                scientific_model_changed=False,
                surface_grid_changed=False,
                mode="fused_diagonal_and_CPU_primitive_surface_blocks",
            )
            atomic_write_json(path, record)
            return original(path, geometry)

        projected.ensure_projected_thermal = seal
        try:
            return original_runner()
        finally:
            projected.ensure_projected_thermal = original

    base.runner = runner
    return base


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=("prepare", "replay", "coordinate", "worker", "launch", "collect"),
    )
    parser.add_argument("--mult", type=int, choices=(1, 2))
    parser.add_argument("--campaign", type=Path, default=CAMPAIGN)
    parser.add_argument("--stage", choices=("energy", "thermal"), default="thermal")
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--max-hours", type=float, default=48)
    args = parser.parse_args()
    if args.campaign.resolve() != CAMPAIGN or args.max_hours <= 0:
        parser.error("Use this frozen campaign and a positive time window")
    module = configure()
    if args.action == "prepare":
        module.prepare()
        plan = api.read_plan(CAMPAIGN)
        plan["producer_sha256"].update(
            {str(MODULE): api.digest(MODULE), str(SOURCE): api.digest(SOURCE)}
        )
        plan["execution_change"] = (
            "Fused diagonal PCM contraction and bounded unchanged CPU curvature "
            "primitive; full replays and physical gates required"
        )
        atomic_write_json(CAMPAIGN / "plan.json", plan)
        api.collect(CAMPAIGN)
    elif args.action == "replay":
        if args.mult is None:
            parser.error("replay requires --mult")
        module.replay(args.mult)
    elif args.action == "coordinate":
        module.coordinate(args.max_hours)
    elif args.action == "collect":
        api.collect(CAMPAIGN)
    else:
        getattr(module.runner(), args.action)(CAMPAIGN, args.stage, args.pilot)
