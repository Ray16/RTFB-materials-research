#!/usr/bin/env python
"""Prepare a self-contained, compiled Overleaf upload without creating a project."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "manuscript"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export_dependencies(text, paper, root, evidence_directory, source_prefix, snapshot, destination, active=()):
    figures = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", text)
    for figure in figures:
        source = (paper/figure).resolve()
        if (not source.is_relative_to(root) or source.suffix != '.pdf'
                or not (figure.startswith(source_prefix) or source.parent == evidence_directory)
                or snapshot['snapshot_sha256'].get(str(source.relative_to(root))) != sha(source)):
            raise ValueError('Figure is not bound by the current evidence snapshot')
        existing = destination/'figures'/source.name
        if existing.exists() and sha(existing) != sha(source):
            raise ValueError('Different figure bytes share an export filename')
        existing.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, existing)
        text = text.replace(figure, 'figures/'+source.name)
    for command, name in re.findall(r"\\(input|include)\{([^}]+)\}", text):
        source = (paper/name).resolve()
        if not source.suffix:
            source = source.with_suffix('.tex')
        if (source.parent != evidence_directory or source.suffix != '.tex'
                or snapshot['snapshot_sha256'].get(str(source.relative_to(root))) != sha(source)):
            raise ValueError('TeX input is not bound by the current evidence snapshot')
        if source in active:
            raise ValueError('Cyclic TeX input dependency')
        child = export_dependencies(source.read_text(), paper, root, evidence_directory,
                                    source_prefix, snapshot, destination, (*active, source))
        target = destination/'inputs'/source.name
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.read_text() != child:
            raise ValueError('Different TeX inputs share an export filename')
        target.write_text(child)
        exported = 'inputs/'+(source.stem if command == 'include' else source.name)
        text = text.replace('\\'+command+'{'+name+'}', '\\'+command+'{'+exported+'}')
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=ROOT.parent / "0_redox_overleaf")
    parser.add_argument("--zip", type=Path, default=PAPER / "overleaf_upload.zip")
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    destination = args.destination.resolve()
    if destination.exists() or args.zip.exists():
        raise RuntimeError("Preserve the existing export. Once connected, pull/merge Overleaf before editing.")
    built = json.loads((PAPER / "build_check.json").read_text())
    if sha(PAPER / "manuscript.pdf") != built["pdf_sha256"]:
        raise ValueError("Draft PDF differs from its checked build")
    if any(sha(ROOT / path) != expected for path, expected in built["source_sha256"].items()):
        raise ValueError("Draft source changed after its checked build")
    pointer = json.loads((PAPER / "evidence/current.json").read_text())
    snapshot = json.loads((ROOT / pointer["manifest"]).read_text())
    if any(sha(ROOT / path) != expected for path, expected in snapshot["snapshot_sha256"].items()):
        raise ValueError("Frozen draft evidence changed")
    destination.mkdir()
    (destination / "figures").mkdir()
    source_prefix = "../" + pointer["figures"] + "/"
    evidence_directory = (ROOT / pointer["manifest"]).parent.resolve()
    document_names = {"manuscript": "main", "supporting_information": "supporting_information"}
    for source_name, exported_name in document_names.items():
        if sha(PAPER / (source_name + ".pdf")) != built["documents"][source_name]["pdf_sha256"]:
            raise ValueError("Document PDF differs from checked build: " + source_name)
        text = (PAPER / (source_name + ".tex")).read_text()
        text = export_dependencies(text, PAPER, ROOT, evidence_directory, source_prefix, snapshot, destination)
        (destination / (exported_name + ".tex")).write_text(text)
    shutil.copyfile(PAPER / "references.bib", destination / "references.bib")
    shutil.copyfile(PAPER / "novelty_references.bib", destination / "novelty_references.bib")
    (destination / ".gitignore").write_text("*.aux\n*.bbl\n*.blg\n*.log\n*.out\n*.toc\n*.fls\n*.fdb_latexmk\n*.synctex.gz\nmain.pdf\n")
    (destination / "README.md").write_text(
        "# Redox polymer manuscript\n\nMain document: `main.tex`. Supporting Information: `supporting_information.tex`. Compiler: pdfLaTeX.\n\n"
        "This is a working draft with explicit unfinished analyses. Figures, bibliography and frozen numerical TeX inputs are included.\n\n"
        "The source and frozen computational evidence are maintained in the 0_redox project. "
        "Do not replace collaborators' edits with a fresh export: fetch and merge before each edit and push.\n")
    outputs = {str(p.relative_to(destination)): sha(p) for p in destination.rglob("*") if p.is_file()}
    build_root = Path(os.environ.get("REDOX_OVERLEAF_BUILD_ROOT",
                                    ROOT / "calcs/validation/overleaf_builds")).expanduser().resolve()
    build_root.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix="redox_overleaf_build_", dir=build_root))
    environment = dict(os.environ, BIBINPUTS=str(destination) + os.pathsep + os.environ.get("BIBINPUTS", ""))
    from pypdf import PdfReader
    exported_documents = {}
    for source_name, exported_name in document_names.items():
        command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-output-directory=" + str(temporary), exported_name + ".tex"]
        for index in range(3):
            with (temporary / f"{exported_name}_pass{index}.log").open("w") as log:
                subprocess.run(command, cwd=destination, env=environment, stdout=log, stderr=subprocess.STDOUT, check=True)
            if index == 0:
                with (temporary / (exported_name + "_bibtex.log")).open("w") as log:
                    subprocess.run(["bibtex", exported_name], cwd=temporary, env=environment,
                                   stdout=log, stderr=subprocess.STDOUT, check=True)
        reference = PdfReader(PAPER / (source_name + ".pdf"))
        exported = PdfReader(temporary / (exported_name + ".pdf"))
        if [p.extract_text() for p in reference.pages] != [p.extract_text() for p in exported.pages]:
            raise ValueError("Self-contained export differs from checked document: " + source_name)
        exported_documents[exported_name] = len(exported.pages)
    if any(sha(destination / path) != expected for path, expected in outputs.items()):
        raise ValueError("Export changed during compilation")
    # The receiving service requires ZIP. Store already compressed PDFs and
    # compress only the small text members with one CPU thread.
    with zipfile.ZipFile(args.zip, "x") as archive:
        for name in sorted(outputs):
            compression = zipfile.ZIP_STORED if name.endswith(".pdf") else zipfile.ZIP_DEFLATED
            archive.write(destination / name, name, compress_type=compression, compresslevel=1)
    with zipfile.ZipFile(args.zip) as archive:
        if archive.testzip() is not None or set(archive.namelist()) != set(outputs):
            raise ValueError("Upload archive integrity/count mismatch")
        if any(hashlib.sha256(archive.read(name)).hexdigest() != expected for name, expected in outputs.items()):
            raise ValueError("Upload archive member differs from source")
    proof = dict(status="compiled_self_contained_upload_ready_no_remote_connection", project_created=False,
                 checked_UTC=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), folder=str(destination),
                 archive=str(args.zip.resolve()), archive_sha256=sha(args.zip), members=outputs, pages=exported_documents["main"], documents=exported_documents,
                 compilation_text_matches_checked_draft=True, build_logs=str(temporary),
                 source_sha256=sha(PAPER / "manuscript.tex"), exporter_sha256=sha(__file__))
    (PAPER / "overleaf_export.json").write_text(json.dumps(proof, indent=2) + "\n")
    print(json.dumps(dict(status=proof["status"], zip=str(args.zip), files=len(outputs), pages=exported_documents["main"], documents=exported_documents)))


if __name__ == "__main__":
    main()
