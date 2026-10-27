#!/usr/bin/env python3
"""Build the scientific story deck and bind its actual inputs and rendered PDF."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

from pypdf import PdfReader, PdfWriter

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--review', type=Path, required=True)
    a = p.parse_args()
    work, review = a.work.resolve(), a.review.resolve()
    work.mkdir(parents=True, exist_ok=False)
    review.mkdir(parents=True, exist_ok=True)
    figure_manifest = ROOT/'results/figures/presentation/story/manifest.json'
    figures = json.loads(figure_manifest.read_text())
    evidence = ROOT/figures['evidence']
    if sha(evidence) != figures['evidence_sha256']:
        raise ValueError('Presentation evidence changed')
    plotter = ROOT/'scripts/plotting/presentation/plot_redox_story.py'
    if sha(plotter) != figures['producer_sha256']:
        raise ValueError('Figure producer changed')
    for path, digest in figures['outputs'].items():
        if sha(ROOT/path) != digest:
            raise ValueError('Figure changed: '+path)
    source_paths = [ROOT/'slides/redox_workflow.tex',
                    *sorted((ROOT/'slides/story').glob('*.tex')),
                    Path(__file__).resolve(), figure_manifest, plotter, evidence,
                    ROOT/'results/reviews/2026-10-26/local_motif_complete/summary.json',
                    ROOT/'results/validation/empty_chain_twenty_2026-10-27/comparison.json']
    frozen = {str(path.relative_to(ROOT)): sha(path) for path in source_paths}
    for i in (1, 2):
        with (work/f'pass{i}.log').open('w') as stream:
            subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
                            '-file-line-error', '-recorder', '-output-directory', str(work),
                            'redox_workflow.tex'], cwd=ROOT/'slides', stdout=stream,
                           stderr=subprocess.STDOUT, check=True, timeout=120)
    log = (work/'pass2.log').read_text()
    errors = [line for line in log.splitlines() if any(token in line for token in
              ('Overfull', 'Undefined control sequence', 'LaTeX Error', 'Missing character:'))]
    if errors:
        raise ValueError(errors)
    pdf = work/'redox_workflow.pdf'
    reader = PdfReader(pdf)
    texts = [page.extract_text() or '' for page in reader.pages]
    if len(texts) != 17 or not all(term in texts[3] for term in ('UMA', 'DFT', 'Grafting')):
        raise ValueError('Expected 13 main / 4 supporting slides and explicit UMA pipeline')
    if not all(term in texts[7] for term in ('Grafting', 'DBU', '57.6', '45.1')):
        raise ValueError('Grafting result missing from slide 8')
    for path, digest in frozen.items():
        if sha(ROOT/path) != digest:
            raise ValueError('Build input changed: '+path)
    for line in (work/'redox_workflow.fls').read_text().splitlines():
        if line.startswith('INPUT '):
            path = Path(line[6:])
            if not path.is_absolute():
                path = ROOT/'slides'/path
            path = path.resolve()
            if path.is_file() and path.is_relative_to(ROOT) and path.suffix in ('.tex', '.pdf', '.png'):
                frozen[str(path.relative_to(ROOT))] = sha(path)
    for name in ('redox_workflow.pdf', 'colleague_workflow.pdf'):
        shutil.copyfile(pdf, ROOT/'slides'/name)
    writer = PdfWriter()
    writer.add_page(reader.pages[3])
    pipeline = ROOT/'results/figures/presentation/story/uma_pipeline.pdf'
    with pipeline.open('wb') as stream:
        writer.write(stream)
    record = dict(status='compiler_and_content_checked_visual_review_separate', pages=17,
                  main_slides=13, backup_slides=4, pdf_sha256=sha(pdf),
                  source_sha256=frozen, compiler_errors=errors, work=str(work),
                  pipeline_figure=str(pipeline.relative_to(ROOT)), pipeline_sha256=sha(pipeline),
                  pipeline_slide=4, grafting_deltaG_slide=8, quantum_jobs_submitted=0)
    (review/'build_check.json').write_text(json.dumps(record, indent=2)+'\n')
    (review/'slide_text.json').write_text(json.dumps(texts, indent=2)+'\n')
    print(json.dumps({k: v for k, v in record.items() if k != 'source_sha256'}))


if __name__ == '__main__':
    main()
