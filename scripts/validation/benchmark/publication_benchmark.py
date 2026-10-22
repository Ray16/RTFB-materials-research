#!/usr/bin/env python
"""Freeze a sourced experimental panel or report an independently admitted review."""

import argparse
import json
from pathlib import Path

from redox.core.common import ROOT, load_config
from redox.core.io import atomic_write_json
from redox.validation.benchmark_panel import freeze, review_report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "report"))
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--review", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.action == "freeze":
        if args.panel.exists():
            parser.error("Frozen panel already exists; use a new version")
        sources = [
            ROOT / "config/benchmark_mecn.py",
            *sorted((ROOT / "data/raw/validation/two_wave_mecn").glob("*.csv")),
        ]
        panel = freeze(load_config("benchmark_mecn").BENCHMARK, sources)
        atomic_write_json(args.panel, panel)
        print(
            json.dumps(
                dict(
                    panel=str(args.panel),
                    targets=len(panel["targets"]),
                    molecules=len({r["id"] for r in panel["targets"]}),
                )
            )
        )
    else:
        if args.review is None or args.out is None:
            parser.error("report requires --review and --out")
        report = review_report(json.loads(args.panel.read_text()), args.review)
        atomic_write_json(args.out, report)
        print(
            json.dumps(
                {
                    k: report[k]
                    for k in ("status", "target_events", "admitted_events", "metrics")
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
