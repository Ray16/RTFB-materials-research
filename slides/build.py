#!/usr/bin/env python3
"""Build the focused deck with checked source data and durable compiler records."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


ROOT=Path(__file__).resolve().parents[1]
SLIDES=ROOT/'slides'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    pointer=json.loads((ROOT/'manuscript/evidence/current.json').read_text());epoch=(ROOT/pointer['manifest']).parent
    spec=importlib.util.spec_from_file_location('deck_checked_directional_values',ROOT/'manuscript/check_tempo_solution.py')
    checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker);descriptor=checker.audit(epoch)
    spec=importlib.util.spec_from_file_location('deck_checked_viologen_values',ROOT/'manuscript/check_viologen_native.py')
    viologen=importlib.util.module_from_spec(spec);spec.loader.exec_module(viologen)
    redox=viologen.audit(epoch,ROOT/pointer['figures'])
    spec=importlib.util.spec_from_file_location('deck_checked_viologen_graft',ROOT/'manuscript/check_viologen_grafting.py')
    grafting=importlib.util.module_from_spec(spec);spec.loader.exec_module(grafting);graft=grafting.audit(epoch)
    paths=[SLIDES/'redox_workflow.tex',SLIDES/'story/motivation_results.tex',SLIDES/'story/supporting.tex',
           epoch/'tempo_solution_values.tex',ROOT/'manuscript/check_tempo_solution.py',
           epoch/'viologen_native_values.tex',ROOT/'manuscript/check_viologen_native.py',Path(__file__).resolve()]
    if graft is not None:paths.extend([epoch/'viologen_grafting_values.tex',ROOT/'manuscript/check_viologen_grafting.py'])
    inputs={str(p.relative_to(ROOT)):sha(p) for p in paths}
    for source in paths[:3]:
        text=source.read_text()
        for reference in re.findall(r'\{([^{}]*\.pdf)\}',text):
            path=(SLIDES/reference).resolve()
            if not path.is_relative_to(ROOT/'results/figures') or not path.exists():
                raise ValueError('Explicit local figure missing or outside its collection')
            inputs[str(path.relative_to(ROOT))]=sha(path)
    build_root=ROOT/'calcs/validation/slides_builds';build_root.mkdir(parents=True,exist_ok=True)
    build=Path(tempfile.mkdtemp(prefix='redox_story_',dir=build_root))
    argv=['pdflatex','-interaction=nonstopmode','-halt-on-error','-output-directory='+str(build),'redox_workflow.tex']
    for index in range(2):
        with (build/f'pass{index+1}.log').open('w') as log:subprocess.run(argv,cwd=SLIDES,stdout=log,stderr=subprocess.STDOUT,check=True)
    log=(build/'redox_workflow.log').read_text()
    if 'Undefined control sequence' in log or 'There were undefined' in log:raise ValueError('Unresolved deck input')
    from pypdf import PdfReader
    pdf=build/'redox_workflow.pdf';reader=PdfReader(pdf);main_frames=(SLIDES/'story/motivation_results.tex').read_text().count('\\begin{frame}')
    supporting_frames=(SLIDES/'story/supporting.tex').read_text().count('\\begin{frame}')
    if len(reader.pages)!=main_frames+supporting_frames:raise ValueError('Deck frame inventory differs')
    text='\n'.join(page.extract_text() for page in reader.pages)
    for number in (f"{descriptor['reduction_eV']:.3f}",f"{descriptor['oxidation_eV']:.3f}"):
        if number not in text:raise ValueError('Checked directional value absent from the deck')
    for value in redox['potentials_V_vs_Fc']:
        if f'{abs(value):.3f}' not in text:raise ValueError('Checked viologen potential absent from the deck')
    if graft is not None and f"{abs(graft['deltaG_kJ_mol_by_cutoff']['50']):.2f}" not in text:
        raise ValueError('Checked exact-viologen grafting value absent from the deck')
    if any(sha(ROOT/path)!=digest for path,digest in inputs.items()):raise ValueError('Deck source changed during compilation')
    shutil.copyfile(pdf,SLIDES/'redox_workflow.pdf')
    report=dict(status='focused_deck_built_checked_directional_result',pages=len(reader.pages),main_frames=main_frames,
        supporting_frames=supporting_frames,source_sha256=inputs,pdf_sha256=sha(SLIDES/'redox_workflow.pdf'),
        build_logs=str(build),overfull_boxes=[line for line in log.splitlines() if line.startswith('Overfull ')],
        visual_review='pending; record separately bound to this PDF',submission_ready=False)
    (SLIDES/'build_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
