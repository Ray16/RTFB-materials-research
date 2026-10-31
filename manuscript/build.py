#!/usr/bin/env python
"""Build the working paper only from unchanged, frozen draft evidence."""

import hashlib
import importlib.util
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


def compiled_figure_inputs(paper,root,snapshot):
    inputs={}
    for document in ('manuscript.tex','supporting_information.tex'):
        for name in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',(paper/document).read_text()):
            figure=(paper/name).resolve();relative=str(figure.relative_to(root))
            if snapshot.get(relative)!=sha(figure):
                raise ValueError('Compiled figure is not bound by the current evidence: '+relative)
            inputs[relative]=sha(figure)
    return inputs


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
    figure_inputs=compiled_figure_inputs(PAPER,ROOT,evidence['snapshot_sha256'])
    from check_claims import audit

    claims = audit()
    ionic_preparation = manifest.parent/'ionic_preparation_replay.json'
    if ionic_preparation.exists():
        expected=json.loads(ionic_preparation.read_text())
        checker_path=ROOT/'scripts/validation/check_charged_polymer_preparation.py'
        definition=importlib.util.spec_from_file_location('checked_ionic_preparation',checker_path)
        checker=importlib.util.module_from_spec(definition);definition.loader.exec_module(checker)
        if checker.check(Path(expected['source']))!=expected:
            raise ValueError('Frozen ionic initial-structure audit differs from actual evidence')
    # Recent context references are separate from the frozen 48-reference base
    # audit. Their primary DOI records must match the actual compiled BibTeX.
    citation_source = ROOT/'scripts/validation/check_recent_polymer_literature.py'
    citation_registry = ROOT/'validation/route_evidence_followup_2026-10-30/recent_references.json'
    definition = importlib.util.spec_from_file_location('checked_recent_polymer_citations',citation_source)
    citations = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(citations)
    citation_check = citations.audit(PAPER/'novelty_references.bib',citation_registry)
    solution_citation_registry = manifest.parent/'tempo_solution_method_references.json'
    solution_citation_check = citations.audit(PAPER/'novelty_references.bib',solution_citation_registry) if solution_citation_registry.exists() else None
    for program in ("pdflatex", "bibtex"):
        if shutil.which(program) is None:
            raise RuntimeError(f"Missing manuscript build program: {program}")
    source = {str(p.relative_to(ROOT)): sha(p) for p in
              (PAPER / "manuscript.tex", PAPER / "supporting_information.tex", PAPER / "references.bib", PAPER / "novelty_references.bib", PAPER / "check_claims.py", PAPER / "check_polymer_replay.py", PAPER / "check_polymer_continuations.py", PAPER / "check_empty_batches.py", PAPER / "check_native_cache.py", PAPER / "check_tempo_solution.py", PAPER / "check_viologen_native.py", Path(__file__), manifest, pointer)}
    source.update(figure_inputs)
    # Bind all frozen inputs through the end of compilation, including cached
    # raw structure records used by the prospective ionic-model demonstration.
    source.update(evidence['snapshot_sha256'])
    if (manifest.parent/'viologen_grafting_result.json').exists():
        source[str((PAPER/'check_viologen_grafting.py').relative_to(ROOT))]=sha(PAPER/'check_viologen_grafting.py')
    if (manifest.parent/'exact_molecular_CPU_validation.json').exists():
        example=ROOT/'examples/exact_viologen_replay';package=example/'data'
        for path in [PAPER/'check_molecular_example.py',example/'replay.py',package/'manifest.json',
                     *[package/name for name in json.loads((package/'manifest.json').read_text())['files_sha256']]]:
            source[str(path.relative_to(ROOT))]=sha(path)
    for path in (citation_source,citation_registry,*[citation_registry.parent/r['metadata_file']
                 for r in json.loads(citation_registry.read_text())['records']]):
        source[str(path.relative_to(ROOT))]=sha(path)
    if solution_citation_check:
        source[str(solution_citation_registry.relative_to(ROOT))]=sha(solution_citation_registry)
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
                  additional_citation_check=citation_check,
                  solution_model_citation_check=solution_citation_check,
                  visual_review="pending; record separately in visual_review.json, bound to these PDF hashes")
    (PAPER / "claim_check.json").write_text(json.dumps(claims, indent=2) + "\n")
    (PAPER / "build_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
