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
    editorial_figure_check=None
    if (manifest.parent/'editorial_figure_review.json').exists():
        from check_editorial_figures import audit as audit_editorial
        editorial_figure_check=audit_editorial(manifest.parent,ROOT/json.loads(pointer.read_text())['figures'])
    screening_check = None
    if (manifest.parent/'screening_V100_integration.json').exists():
        from check_screening_controls import audit as audit_screening
        screening_check = audit_screening(manifest.parent)
    ionic_solver_check = None
    if (manifest.parent/'ionic_compact_review.json').exists():
        from check_charged_pilot import audit as audit_ionic_solver
        ionic_solver_check = audit_ionic_solver(manifest.parent)
    xtb_runtime_check = None
    if (manifest.parent/'xtb_runtime_installation.json').exists():
        from check_xtb_runtime import audit as audit_xtb_runtime
        xtb_runtime_check = audit_xtb_runtime(manifest.parent)
    context_consequences_check = None
    if (manifest.parent/'context_consequences.json').exists():
        from check_context_consequences import audit as audit_context_consequences
        context_consequences_check = audit_context_consequences(manifest.parent)
    component_sensitivity_check = None
    if (manifest.parent/'component_implementation_inputs.json').exists():
        from check_component_sensitivity import audit as audit_component_sensitivity
        component_sensitivity_check = audit_component_sensitivity(manifest.parent)
    native_component_check = None
    if (manifest.parent/'ionic_native_component_review.json').exists():
        from check_native_component_pilot import audit as audit_native_component
        native_component_check = audit_native_component(manifest.parent)
    ionic_native_cycle_check = None
    if (manifest.parent/'ionic_native_cycles.json').exists():
        from check_ionic_native_cycles import audit as audit_ionic_native_cycles
        ionic_native_cycle_check = audit_ionic_native_cycles(manifest.parent)
    motif_interaction_check = None
    if (manifest.parent/'motif_energy_collection.json').exists():
        from check_motif_interactions import audit as audit_motif_interactions
        motif_interaction_check = audit_motif_interactions(manifest.parent)
    counterion_association_check = None
    if (manifest.parent/'counterion_inputs.json').exists():
        from check_counterion_association import audit as audit_counterion_association
        counterion_association_check = audit_counterion_association(manifest.parent)
    registry_expansion_check = None
    if (manifest.parent/'registry_update_current.json').exists():
        from check_registry_expansion import audit as audit_registry_expansion
        registry_expansion_check = audit_registry_expansion(manifest.parent)
    occupied_sampling_check = None
    if (manifest.parent/'occupied_sampling_manifest.json').exists():
        from redox.validation.sampling_publication import audit as audit_occupied_sampling
        occupied_sampling_check = audit_occupied_sampling(manifest.parent, prefix='occupied_sampling_')
    primary_native_cohort_check = None
    corrected_example_check = None
    if (manifest.parent/'primary_native_cohort_inputs.json').exists():
        from check_primary_native_cohort import audit as audit_primary_native_cohort
        primary_native_cohort_check = audit_primary_native_cohort(manifest.parent)
        from check_corrected_example import audit as audit_corrected_example
        corrected_example_check = audit_corrected_example(ROOT,manifest.parent)
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
    registry_citation_registry = manifest.parent/'registry_update_reference_registry.json'
    registry_citation_check = citations.audit(PAPER/'novelty_references.bib',registry_citation_registry) if registry_citation_registry.exists() else None
    for program in ("pdflatex", "bibtex"):
        if shutil.which(program) is None:
            raise RuntimeError(f"Missing manuscript build program: {program}")
    source = {str(p.relative_to(ROOT)): sha(p) for p in
              (PAPER / "manuscript.tex", PAPER / "supporting_information.tex", PAPER / "references.bib", PAPER / "novelty_references.bib", PAPER / "check_claims.py", PAPER / "check_polymer_replay.py", PAPER / "check_polymer_continuations.py", PAPER / "check_empty_batches.py", PAPER / "check_native_cache.py", PAPER / "check_tempo_solution.py", PAPER / "check_viologen_native.py", Path(__file__), manifest, pointer)}
    source.update(figure_inputs)
    if editorial_figure_check is not None:source['manuscript/check_editorial_figures.py']=sha(PAPER/'check_editorial_figures.py')
    if native_component_check is not None:
        source['manuscript/check_native_component_pilot.py']=sha(PAPER/'check_native_component_pilot.py')
    if ionic_native_cycle_check is not None:
        source['manuscript/check_ionic_native_cycles.py']=sha(PAPER/'check_ionic_native_cycles.py')
        source['src/redox/workflows/ionic_native_cycles.py']=sha(ROOT/'src/redox/workflows/ionic_native_cycles.py')
    if motif_interaction_check is not None:
        source['manuscript/check_motif_interactions.py']=sha(PAPER/'check_motif_interactions.py')
    if counterion_association_check is not None:
        for name in ('manuscript/check_counterion_association.py','src/redox/properties/counterion.py',
                     'src/redox/properties/context.py','src/redox/properties/conformer_geometry.py',
                     'src/redox/properties/conformer_ensemble.py'):
            source[name]=sha(ROOT/name)
    if registry_expansion_check is not None:
        for name in ('manuscript/check_registry_expansion.py','scripts/pipeline/screen_grafting_routes.py',
                     'src/redox/build/attachment_identity.py'):
            source[name]=sha(ROOT/name)
    if occupied_sampling_check is not None:
        for name in ('sampling_publication.py','sampling_batches.py','native_conformer_batches.py'):
            relative='src/redox/validation/'+name; source[relative]=sha(ROOT/relative)
        for relative in ('src/redox/properties/conformer_geometry.py','src/redox/properties/conformer_ensemble.py',
                         'scripts/plotting/plot_occupied_sampling.py','results/figures/slides/occupied_sampling_convergence.pdf'):
            source[relative]=sha(ROOT/relative)
    if primary_native_cohort_check is not None:
        source['manuscript/check_primary_native_cohort.py']=sha(PAPER/'check_primary_native_cohort.py')
        source['src/redox/validation/native_component_sensitivity.py']=sha(ROOT/'src/redox/validation/native_component_sensitivity.py')
        example=ROOT/'examples/exact_viologen_corrected_replay';cache=example/'data'
        for path in [PAPER/'check_corrected_example.py',example/'replay.py',cache/'manifest.json',
                     *[cache/name for name in json.loads((cache/'manifest.json').read_text())['files_sha256']]]:
            source[str(path.relative_to(ROOT))]=sha(path)
    if component_sensitivity_check is not None:
        source['manuscript/check_component_sensitivity.py']=sha(PAPER/'check_component_sensitivity.py')
        source['src/redox/validation/component_sensitivity.py']=sha(ROOT/'src/redox/validation/component_sensitivity.py')
    if context_consequences_check is not None:
        source['manuscript/check_context_consequences.py']=sha(PAPER/'check_context_consequences.py')
        source['src/redox/validation/context_consequences.py']=sha(ROOT/'src/redox/validation/context_consequences.py')
    if screening_check is not None:
        source['manuscript/check_screening_controls.py']=sha(PAPER/'check_screening_controls.py')
        if 'electronic_benchmark_implementation_sensitivity' in screening_check:
            source['src/redox/validation/screening_sensitivity.py']=sha(ROOT/'src/redox/validation/screening_sensitivity.py')
    if ionic_solver_check is not None:source['manuscript/check_charged_pilot.py']=sha(PAPER/'check_charged_pilot.py')
    if xtb_runtime_check is not None:source['manuscript/check_xtb_runtime.py']=sha(PAPER/'check_xtb_runtime.py')
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
        if occupied_sampling_check is not None:
            sampling_text='\n'.join(page.extract_text() for page in reader.pages)
            numerical=re.sub(r'(?<=\d)\s*\.\s*(?=\d)', '.', sampling_text)
            required=['2.40','3.88'] if name=='manuscript' else ['2.55','2.40','2.54','0.40','0.65','1.28','3.88']
            if any(value not in numerical for value in required):
                raise ValueError('Checked occupied sampling values absent from '+name)
        if ionic_native_cycle_check is not None:
            from check_ionic_native_cycles import check_compiled_values
            check_compiled_values('\n'.join(page.extract_text() for page in reader.pages),ionic_native_cycle_check)
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
                  screening_implementation_check=screening_check,
                  ionic_solver_check=ionic_solver_check,
                  context_consequences_check=context_consequences_check,
                  component_sensitivity_check=component_sensitivity_check,
                  native_component_check=native_component_check,
                  ionic_native_cycle_check=ionic_native_cycle_check,
                  motif_interaction_check=motif_interaction_check,
                  counterion_association_check=counterion_association_check,
                  registry_expansion_check=registry_expansion_check,
                  editorial_figure_check=editorial_figure_check,
                  occupied_sampling_check=occupied_sampling_check,
                  registry_citation_check=registry_citation_check,
                  primary_native_cohort_check=primary_native_cohort_check,
                  corrected_example_check=corrected_example_check,
                  xtb_runtime_check=xtb_runtime_check,
                  visual_review="pending; record separately in visual_review.json, bound to these PDF hashes")
    (PAPER / "claim_check.json").write_text(json.dumps(claims, indent=2) + "\n")
    (PAPER / "build_check.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
