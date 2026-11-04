#!/usr/bin/env python3
"""Build the Overleaf article, SI and TOC with capped local CPU use."""
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from pypdf import PdfReader
from check_narrative import audit

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', BIBINPUTS=str(ROOT) + os.pathsep)
    checks = ROOT.parents[1] / 'calcs/validation/overleaf_checks'
    checks.mkdir(exist_ok=True)
    build = Path(tempfile.mkdtemp(prefix='narrative_revision_', dir=checks))
    toc = ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory=' + str(build), 'graphical_abstract.tex']
    with (build / 'toc_build.log').open('w') as log:
        subprocess.run(toc, cwd=ROOT / 'figures', env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    toc_pdf = PdfReader(build / 'graphical_abstract.pdf')
    assert len(toc_pdf.pages) == 1, 'TOC overflows one page'
    shutil.copyfile(build / 'graphical_abstract.pdf', ROOT / 'figures/graphical_abstract.pdf')
    claim_check = audit()
    (ROOT / 'narrative_claim_check.json').write_text(json.dumps(claim_check, indent=2) + '\n')
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in ROOT.rglob('*') if p.is_file() and '.git' not in p.parts and p.suffix in {'.tex', '.bib', '.csv', '.pdf'} and p.name not in {'main.pdf', 'supporting_information.pdf'}}
    documents = {}
    for name in ['main', 'supporting_information']:
        command = ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory=' + str(build), name + '.tex']
        for i in range(6):
            with (build / f'{name}_pass{i}.log').open('w') as log:
                subprocess.run(command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            if i == 0:
                with (build / f'{name}_bibtex.log').open('w') as log:
                    subprocess.run(['bibtex', name], cwd=build, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            if i >= 2:
                current_log = (build / (name + '.log')).read_text(errors='replace')
                if 'undefined' not in current_log.lower() and 'Rerun to get cross-references right' not in current_log:
                    break
        log_text = (build / (name + '.log')).read_text(errors='replace')
        unresolved = [line for line in log_text.splitlines() if 'undefined' in line.lower() or 'Rerun to get cross-references right' in line]
        assert not unresolved, (name, unresolved)
        overfull = [line for line in log_text.splitlines() if 'Overfull' in line]
        reader = PdfReader(build / (name + '.pdf'))
        text = '\n'.join(page.extract_text() for page in reader.pages)
        assert 'Working draft' in text
        assert all('�' not in page.extract_text() for page in reader.pages)
        documents[name] = {'pages': len(reader.pages), 'pdf_sha256': sha(build / (name + '.pdf')), 'overfull_boxes': overfull, 'unresolved_references': unresolved}
        print(name, documents[name], flush=True)
    assert all(sha(ROOT / name) == digest for name, digest in inputs.items()), 'Document input changed during build'
    result = {'status': 'compiled_with_reconstructed_benchmark_and_unchanged_frozen_inputs', 'checked_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'build_directory': str(build), 'documents': documents, 'input_sha256': inputs, 'claim_check_sha256': sha(ROOT / 'narrative_claim_check.json'), 'new_quantum_calculations': False, 'submission_ready': False, 'visual_review': 'pending'}
    (ROOT / 'narrative_build_check.json').write_text(json.dumps(result, indent=2) + '\n')
    print('build_directory', build, flush=True)

if __name__ == '__main__':
    main()
