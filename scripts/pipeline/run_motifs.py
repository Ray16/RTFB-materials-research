#!/usr/bin/env python
"""Compatibility command; maintained orchestration lives in redox.workflows.motifs."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from redox.workflows.motifs import *  # noqa: F403
from redox.workflows import motifs


def resume_config_matches(previous, config, scaffold, groups):
    return motifs.resume_config_matches(previous, config, scaffold, groups, root=ROOT)


if __name__ == "__main__":
    main()
