#!/usr/bin/env python
"""Run the unchanged full spin control with the XC temporary-lifetime adapter."""
import importlib.util
import hashlib
import json
from pathlib import Path
import sys


def main():
    root = Path(__file__).resolve().parents[2]
    original = root / "scripts/validation/qualify_smd_workspace.py"
    spec = importlib.util.spec_from_file_location("XC_full_spin_control", original)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    before = module.adapted_body

    def adapted(body):
        text = before(body)
        old = "scripts/pipeline/smd_workspace_execution.py"
        if text.count(old) != 1:
            raise ValueError("Full workspace-control template changed")
        return text.replace(old, "scripts/pipeline/xc_workspace_execution.py")

    module.adapted_body = adapted
    try:
        module.main()
    finally:
        path = Path(sys.argv[sys.argv.index("--out") + 1]) / "status.json"
        if path.exists():
            result = json.loads(path.read_text())
            result.setdefault("source_sha256", {}).update({str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (Path(__file__).resolve(), original)})
            result["execution_variant"] = "Release completed RKS density-loop AO batch before XC grid-response work"
            path.write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
