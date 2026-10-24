#!/usr/bin/env python
"""Full original-versus-blocked spin controls for PCM first-derivative storage."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys


def main():
    root = Path(__file__).resolve().parents[2]
    original = root / "scripts/validation/qualify_smd_workspace.py"
    spec = importlib.util.spec_from_file_location("PCM_first_response_spin_control", original)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    before = module.adapted_body

    def adapted(body):
        text = before(body)
        old = "scripts/pipeline/smd_workspace_execution.py"
        if text.count(old) != 1:
            raise ValueError("Full spin-control adapter template changed")
        return text.replace(old, "scripts/pipeline/smd_first_response_workspace.py")

    module.adapted_body = adapted
    try:
        module.main()
    finally:
        path = Path(sys.argv[sys.argv.index("--out")+1]) / "status.json"
        if path.exists():
            result = json.loads(path.read_text())
            result.setdefault("source_sha256", {}).update({str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (Path(__file__).resolve(), original,
                    root / "scripts/pipeline/smd_first_response_workspace.py",
                    root / "scripts/pipeline/pcm_first_derivative_blocks.py")})
            result["execution_variant"] = "Bounded first-derivative rectangles from the unchanged installed pairwise PCM primitive"
            path.write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()
