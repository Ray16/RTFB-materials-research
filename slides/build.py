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
    spec=importlib.util.spec_from_file_location('deck_checked_ionic_native',ROOT/'manuscript/check_native_component_pilot.py')
    ionic=importlib.util.module_from_spec(spec);spec.loader.exec_module(ionic);native=ionic.audit(epoch)
    spec=importlib.util.spec_from_file_location('deck_checked_primary_native_cohort',ROOT/'manuscript/check_primary_native_cohort.py')
    primary=importlib.util.module_from_spec(spec);spec.loader.exec_module(primary);cohort=primary.audit(epoch)
    ionic_cycle=None
    if (epoch/'ionic_native_cycles.json').exists():
        spec=importlib.util.spec_from_file_location('deck_checked_ionic_cycle',ROOT/'manuscript/check_ionic_native_cycles.py')
        cycle_checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(cycle_checker);ionic_cycle=cycle_checker.audit(epoch)
    motif_check=None
    if (epoch/'motif_energy_collection.json').exists():
        spec=importlib.util.spec_from_file_location('deck_checked_motif_interactions',ROOT/'manuscript/check_motif_interactions.py')
        checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker);motif_check=checker.audit(epoch)
    counterion_check=None
    if (epoch/'counterion_inputs.json').exists():
        spec=importlib.util.spec_from_file_location('deck_checked_counterion_association',ROOT/'manuscript/check_counterion_association.py')
        checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker);counterion_check=checker.audit(epoch)
    registry_check=None
    if (epoch/'registry_update_current.json').exists():
        spec=importlib.util.spec_from_file_location('deck_checked_registry_expansion',ROOT/'manuscript/check_registry_expansion.py')
        checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker);registry_check=checker.audit(epoch)
    sampling_folder=SLIDES/'evidence/2026-11-02_occupied_sampling'
    spec=importlib.util.spec_from_file_location('deck_checked_occupied_sampling',SLIDES/'check_sampling.py')
    checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker);sampling_check=checker.audit(sampling_folder)
    paths=[SLIDES/'redox_workflow.tex',SLIDES/'story/motivation_results.tex',SLIDES/'story/supporting.tex',
           epoch/'tempo_solution_values.tex',ROOT/'manuscript/check_tempo_solution.py',
           epoch/'viologen_native_values.tex',ROOT/'manuscript/check_viologen_native.py',Path(__file__).resolve()]
    if graft is not None:paths.extend([epoch/'viologen_grafting_values.tex',ROOT/'manuscript/check_viologen_grafting.py'])
    paths.extend([epoch/'ionic_native_component_review.json',ROOT/'manuscript/check_native_component_pilot.py'])
    paths.extend([epoch/name for name in ('primary_native_cohort_inputs.json','primary_native_cohort_cycles.json','primary_native_component_reviews.json','primary_native_cohort_table.tex')])
    paths.extend([ROOT/'manuscript/check_primary_native_cohort.py',ROOT/'src/redox/validation/native_component_sensitivity.py'])
    if ionic_cycle is not None:
        suffix=next(s for s in ('_v4','_v3','_v2','') if (epoch/('ionic_native_cycles'+s+'.json')).exists())
        paths.extend([ROOT/'manuscript/check_ionic_native_cycles.py',ROOT/'src/redox/workflows/ionic_native_cycles.py',
            epoch/'ionic_native_cycle_preparation.json',
            *[epoch/(stem+suffix+extension) for stem,extension in [('ionic_native_cycles','.json'),('ionic_native_cycle_values','.tex'),('ionic_native_cycle_table','.tex')]]])
    if motif_check is not None:
        paths.extend([ROOT/'manuscript/check_motif_interactions.py',*[p for p in epoch.glob('motif_*') if p.is_file()]])
    if counterion_check is not None:
        paths.extend([ROOT/'manuscript/check_counterion_association.py',
            *[ROOT/'src/redox/properties'/name for name in ('counterion.py','context.py','conformer_geometry.py','conformer_ensemble.py')],
            *[p for p in epoch.glob('counterion_*') if p.is_file()]])
    if registry_check is not None:
        paths.extend([ROOT/'manuscript/check_registry_expansion.py',ROOT/'scripts/pipeline/screen_grafting_routes.py',
            ROOT/'src/redox/build/attachment_identity.py',*[p for p in epoch.glob('registry_update_*') if p.is_file()]])
    paths.extend([SLIDES/'check_sampling.py',ROOT/'src/redox/validation/sampling_batches.py',
        ROOT/'src/redox/validation/native_conformer_batches.py',ROOT/'src/redox/validation/sampling_publication.py',
        ROOT/'scripts/plotting/plot_occupied_sampling.py',
        *[p for p in sampling_folder.iterdir() if p.is_file()]])
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
    if ionic_cycle is not None:
        cycle_checker.check_compiled_values(text,ionic_cycle)
        count={1:'One',2:'Two',3:'Three',4:'Four'}[ionic_cycle['complete_finite_cycles']]
        if f"{abs(ionic_cycle['delta_delta_G_by_cutoff_kJ_mol']['50']):.2f}" not in text or count+'offour' not in ''.join(text.split()):
            raise ValueError('Checked ionic cycle or its finite coverage is absent from the deck')
    elif (str(native['atoms']) not in text or f"{native['minimum_internal_frequency_cm']:.2f}" not in text
            or 'occupiedfree-energycycle' not in ''.join(text.split())):
        raise ValueError('Checked ionic component or its cycle limit is absent from the deck')
    if motif_check is not None and 'motifrepulsion' not in ''.join(text.split()):
        raise ValueError('Checked motif/whole-chain interpretation is absent from the deck')
    if counterion_check is not None:
        # Math-mode decimal punctuation can acquire extraction-only spaces.
        numerical_text=re.sub(r'(?<=\d)\s*\.\s*(?=\d)', '.',text)
        for number in (f"{counterion_check['minimum_primary_association_1M_kJ_mol']:.2f}",
                       f"{counterion_check['maximum_primary_association_1M_kJ_mol']:.2f}"):
            if number not in numerical_text:raise ValueError('Checked chloride association range absent from the deck')
        if 'threeduplicategroups' not in ''.join(text.lower().split()):
            raise ValueError('Chloride duplicate scope absent from the deck')
    if registry_check is not None and str(registry_check['exact_graphs']) not in text:
        raise ValueError('Checked current registry census absent from the deck')
    sampling_page=reader.pages[18].extract_text()
    if ('failsthesamplingcheck' not in ''.join(sampling_page.lower().split())
            or 'Twenty' not in sampling_page or sampling_check['numbers'][-1] not in sampling_page):
        raise ValueError('Checked occupied sampling outcome absent from supporting slide 19')
    if '81' not in text or cohort['complete_cycles']!=6:raise ValueError('Checked primary-native coverage absent from the deck')
    if any(sha(ROOT/path)!=digest for path,digest in inputs.items()):raise ValueError('Deck source changed during compilation')
    shutil.copyfile(pdf,SLIDES/'redox_workflow.pdf')
    report=dict(status='focused_deck_built_checked_directional_result',pages=len(reader.pages),main_frames=main_frames,
        supporting_frames=supporting_frames,source_sha256=inputs,pdf_sha256=sha(SLIDES/'redox_workflow.pdf'),
        build_logs=str(build),overfull_boxes=[line for line in log.splitlines() if line.startswith('Overfull ')],
        ionic_cycle_check=ionic_cycle,motif_interaction_check=motif_check,counterion_association_check=counterion_check,
        registry_expansion_check=registry_check,sampling_check=sampling_check,
        visual_review='pending; record separately bound to this PDF',submission_ready=False)
    (SLIDES/'build_check.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
