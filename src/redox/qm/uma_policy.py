"""Configuration and explicit admission for the current UMA engines."""
import json
import os
from pathlib import Path

from redox.core.common import ROOT

POLICY = ROOT / "config/uma_execution.json"


def load_policy(path=None):
    path = POLICY if path is None else Path(path)
    policy = json.loads(path.read_text())
    for key, variable in (("nitro_python", "REDOX_NITRO_PYTHON"),
                          ("legacy_python", "REDOX_LEGACY_PYTHON"),
                          ("nitro_library_path", "REDOX_NITRO_LIBRARY_PATH"),
                          ("fairchem_cache", "FAIRCHEM_CACHE_DIR")):
        if os.environ.get(variable):
            policy[key] = os.environ[variable]
    if policy.get("model") != "uma-s-1p2p1" or policy.get("task") != "omol":
        raise ValueError("The current UMA adapter requires UMA-S-1.2.1 / omol")
    for key in ("nitro_output", "legacy_output"):
        output = Path(policy[key])
        if output.is_absolute() or ".." in output.parts:
            raise ValueError(f"{key} must remain inside the current workspace")
    return policy


def require_production_qualification(policy):
    """Do not confuse an installed, executing engine with numerical qualification."""
    certificate = ROOT / policy.get("qualification_file", "config/uma_nitro_qualification.json")
    record = json.loads(certificate.read_text())
    if record.get("production_qualified") is not True:
        raise RuntimeError(
            "UMA Nitro is configured but not production-qualified: retained force "
            "repeatability checks failed. Run redox uma-check and the isolated Nitro "
            "diagnostic; do not reuse diagnostic output as a production warm start."
        )
    # Positive admission requires a source-bound, complete numerical review, not a flag alone.
    from redox.core.hashing import sha256_file
    sources = record.get("source_sha256", {})
    checks = record.get("checks", {})
    required = ("full_panel", "repeatability", "energy_force_parity", "optimizer")
    if (not sources or any(checks.get(key) is not True for key in required)
            or record.get("checkpoint_sha256") != policy["checkpoint_sha256"]
            or record.get("reference_sha256") != policy["reference_sha256"]):
        raise RuntimeError("Incomplete Nitro qualification certificate")
    for path, expected in sources.items():
        if sha256_file(ROOT / path) != expected:
            raise RuntimeError("Nitro qualification source or evidence changed")
