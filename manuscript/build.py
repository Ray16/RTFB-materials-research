#!/usr/bin/env python
"""Build the working paper only from unchanged, frozen draft evidence."""

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "manuscript"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check_draft_labels(text):
    # pdfLaTeX may hyphenate a required word at a physical line boundary.
    joined = re.sub(r"(?<=[A-Za-z])-\s*\n\s*(?=[A-Za-z])", "", text)
    joined = re.sub(r"\s+", " ", joined)
    if (any(label not in joined for label in ("Working draft", "DBU"))
            or re.search(r"\bready for publication\b", joined, re.IGNORECASE)):
        raise RuntimeError("Draft labels, pending analyses or results absent from built PDF")


def main():
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    pointer = PAPER / "evidence/current.json"
    manifest = ROOT / json.loads(pointer.read_text())["manifest"]
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
              (PAPER / "manuscript.tex", PAPER / "supporting_information.tex", PAPER / "references.bib", PAPER / "novelty_references.bib", PAPER / "check_claims.py", PAPER / "check_polymer_replay.py", PAPER / "check_polymer_continuations.py", PAPER / "check_empty_batches.py", Path(__file__), manifest, pointer)}
    # Keep compiler logs on durable project storage rather than node-local /tmp.
    # The root can be configured for an independently relocated checkout.
    build_root = Path(os.environ.get("REDOX_MANUSCRIPT_BUILD_ROOT",
                                    ROOT / "calcs/validation/manuscript_builds")).expanduser().resolve()
    build_root.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix="redox_manuscript_", dir=build_root))
    env = dict(os.environ, BIBINPUTS=str(PAPER) + os.pathsep + os.environ.get("BIBINPUTS", ""))
    documents = {}
    from pypdf import PdfReader
    for name in ("manuscript", "supporting_information"):
        command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-output-directory=" + str(temporary), name + ".tex"]
        for index in range(3):
            with (temporary / f"{name}_pass{index + 1}.log").open("w") as log:
                subprocess.run(command, cwd=PAPER, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            if index == 0:
                with (temporary / (name + "_bibtex.log")).open("w") as log:
                    subprocess.run(["bibtex", name], cwd=temporary, env=env,
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
        log = (temporary / (name + ".log")).read_text()
        if "undefined references" in log or "undefined citations" in log or "There were undefined" in log:
            raise RuntimeError("Unresolved citations or references in " + name)
        pdf = temporary / (name + ".pdf")
        reader = PdfReader(pdf)
        documents[name] = dict(pages=len(reader.pages), pdf_sha256=sha(pdf),
            overfull_boxes=re.findall(r"Overfull \\[hv]box[^\n]*", log))
    if any(sha(ROOT / p) != expected for p, expected in source.items()):
        raise ValueError("Manuscript source changed during build")
    reader = PdfReader(temporary / "manuscript.pdf")
    text = "\n".join(page.extract_text() for page in reader.pages)
    check_draft_labels(text)
    # Publish only after both documents and all source/label checks have passed.
    for name in documents:
        shutil.copyfile(temporary / (name + ".pdf"), PAPER / (name + ".pdf"))
    report = dict(status="draft_built_with_frozen_evidence", checked_UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  pages=len(reader.pages), source_sha256=source, pdf_sha256=sha(PAPER / "manuscript.pdf"),
                  build_logs=str(temporary), submission_ready=False, documents=documents,
                  visual_review="pending; record separately in visual_review.json, bound to these PDF hashes")
    (PAPER / "claim_check.json").write_text(json.dumps(claims, indent=2) + "\n")
    (PAPER / "build_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
