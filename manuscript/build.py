#!/usr/bin/env python
"""Build the working paper only from unchanged, frozen draft evidence."""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "manuscript"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    manifest = PAPER / "evidence/2026-10-21/manifest.json"
    evidence = json.loads(manifest.read_text())
    for path, expected in evidence["snapshot_sha256"].items():
        if sha(ROOT / path) != expected:
            raise ValueError(f"Draft evidence changed: {path}")
    from check_claims import audit

    claims = audit()
    for program in ("pdflatex", "bibtex"):
        if shutil.which(program) is None:
            raise RuntimeError(f"Missing manuscript build program: {program}")
    source = {str(p.relative_to(ROOT)): sha(p) for p in
              (PAPER / "manuscript.tex", PAPER / "references.bib", PAPER / "check_claims.py", Path(__file__), manifest)}
    temporary = Path(tempfile.mkdtemp(prefix="redox_manuscript_"))
    env = dict(os.environ, BIBINPUTS=str(PAPER) + os.pathsep + os.environ.get("BIBINPUTS", ""))
    command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-output-directory=" + str(temporary), "manuscript.tex"]
    for index in range(3):
        with (temporary / f"latex_pass{index + 1}.log").open("w") as log:
            subprocess.run(command, cwd=PAPER, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        if index == 0:
            with (temporary / "bibtex.log").open("w") as log:
                subprocess.run(["bibtex", "manuscript"], cwd=temporary, env=env,
                               stdout=log, stderr=subprocess.STDOUT, check=True)
    log = (temporary / "manuscript.log").read_text()
    if "undefined references" in log or "undefined citations" in log or "There were undefined" in log:
        raise RuntimeError("Unresolved manuscript citations or references")
    if any(sha(ROOT / p) != expected for p, expected in source.items()):
        raise ValueError("Manuscript source changed during build")
    shutil.copyfile(temporary / "manuscript.pdf", PAPER / "manuscript.pdf")
    from pypdf import PdfReader

    reader = PdfReader(PAPER / "manuscript.pdf")
    text = "\n".join(page.extract_text() for page in reader.pages)
    if "Working draft" not in text or "Draft completion item" not in text or "DBU" not in text:
        raise RuntimeError("Draft labels, pending analyses or results absent from built PDF")
    report = dict(status="draft_built_with_frozen_evidence", checked_UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  pages=len(reader.pages), source_sha256=source, pdf_sha256=sha(PAPER / "manuscript.pdf"),
                  build_logs=str(temporary), submission_ready=False, visual_review="pending")
    (PAPER / "claim_check.json").write_text(json.dumps(claims, indent=2) + "\n")
    (PAPER / "build_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
