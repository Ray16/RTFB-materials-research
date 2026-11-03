"""Read-only catalogue of maintained commands, controls and recovery entry points."""
import argparse
import json

CATEGORIES = {
 'production': ['motifs','launch','run-gpu'],
 'design': ['build','build-validation','expand','search','source-audit'],
 'analysis_and_validation': ['status','energy-audit','review','thermal-postprocess','uma-check','seed-reference','resources'],
 'calculation_stages': ['uma','dft','sp','native-thermal'],
 'discovery': ['workflows'],
 'controls': [],
 'recovery': [],
 'historical': [],
}
SPECIALIZED = [
 dict(kind='design', script='scripts/candidates/expand_ring_panel.py',
      purpose='Exact mapped precursor/product ring panel; computes no properties'),
 dict(kind='design', script='scripts/candidates/merge_registry.py',
      purpose='Source-frozen graph/alias union with exact route decisions; preserves existing results'),
 dict(kind='checked_control', module='redox.workflows.sampling_energy',
      purpose='Full 240-slot corrected electronic method-control cohort; separate native and ensemble gates'),
 dict(kind='checked_control', module='redox.workflows.sampling_native',
      purpose='Native method-control inputs require checked exact energies and a completed largest pilot'),
 dict(kind='checked_control', script='scripts/fleet/run_corrected_sampling_queue.py',
      purpose='Maintained resumable electronic coordinator; checked pilots, actual owners and shared resource gate'),
 dict(kind='checked_control', script='scripts/fleet/run_corrected_native_queue.py',
      purpose='Maintained full-roster native coordinator; independent closed review and source-bound followups'),
 dict(kind='checked_analysis', module='redox.validation.sampling_snapshot',
      purpose='Read-only frozen 240-slot receipt/context audit; does not repeat physics or qualify ensembles'),
 dict(kind='checked_control', module='redox.workflows.counterion_geometry',
      purpose='Exact oxidized graft plus chloride (+1 singlet), distinct from the reduced +1 doublet'),
 dict(kind='checked_control', module='redox.workflows.counterion_gradient_pilot',
      purpose='Prospective T4 fixed-gradient comparison only; no optimization, Hessian or native G admission'),
 dict(kind='checked_control', module='redox.workflows.counterion_geometry_resources',
      purpose='T4 geometry-only pilot after checked fixed gradient; full closed review precedes other repairs'),
 dict(kind='checked_control', module='redox.workflows.counterion_thermochemistry',
      purpose='Exact primary electronic/native components and two-to-one solute standard-state conversion'),
 dict(kind='recovery', module='redox.workflows.counterion_mode_repair',
      purpose='Diagnose a failed minimum and prepare both mode signs; fresh geometry/native QC required'),
 dict(kind='recovery', module='redox.workflows.charged_hessian_recovery',
      purpose='Proven closed post-Hessian I/O failure only; reuse verified response/matrix and freshly validate all derivatives'),
 dict(kind='recovery', script='scripts/fleet/run_saved_ionic_validation.py',
      purpose='Bounded saved-ionic-matrix coordinator; exact stopped owner, unchanged physical gates and full resource class'),
 dict(kind='checked_control_preparation', module='redox.workflows.viologen_motif_fragments',
      purpose='Frozen methyl-capped ionic motifs at checked parent frames; no interaction energy from geometry'),
 dict(kind='checked_control', module='redox.workflows.motif_interaction_energy',
      purpose='Matched frozen A/B/AB electronic triplets; no binding G, coupling or population'),
 dict(kind='recovery', module='redox.workflows.charged_partial_recovery',
      purpose='Closed allocation failure only; exact partial frame, full solvent H1 and prospective H100 workspace'),
 dict(kind='checked_analysis', module='redox.workflows.ionic_native_cycles',
      purpose='Read-only matched native occupancy cycles; missing components remain blank'),
 dict(kind='checked_control_preparation', module='redox.workflows.screening_native_panel',
      purpose='Prepare exact primary molecular Hessian repeats; separate pilot launch/review required'),
 dict(kind='checked_campaign', module='redox.workflows.charged_native_followup',
      purpose='Release only source/dimension-matched components after an exact charged-native pilot passes'),
 dict(kind='recovery', module='redox.workflows.relocate_unstarted_namespace',
      purpose='Explicit byte-identical dependency relocation before any scientific start'),
 dict(kind='checked_control', module='redox.workflows.screening_control',
      purpose='Explicit numerical controls; does not change the production method'),
 dict(kind='checked_campaign', module='redox.workflows.charged_polymer_pilot',
      purpose='Route-bound ionic polymer pilots; exact namespace and resource qualification required'),
 dict(kind='checked_campaign', module='redox.workflows.solution_pair_namespace',
      purpose='Matched directional solution response for an exact candidate'),
 dict(kind='recovery', module='redox.workflows.molecular_native_namespace',
      purpose='Private source-bound native continuation; inspect its subcommands and source contract'),
 dict(kind='historical', module='redox.current',
      purpose='Compatibility alias to the public command registry; not a separate implementation'),
]

SPECIALIZED_CATEGORIES = dict(design='design', recovery='recovery', checked_analysis='analysis_and_validation',
    checked_control_preparation='controls', checked_control='controls', checked_campaign='production', historical='historical')


def catalogue(category=None):
    from redox.cli import COMMANDS
    listed = [name for group in CATEGORIES.values() for name in group if name in COMMANDS]
    if len(set(listed)) != len(listed) or set(listed) != set(COMMANDS):
        raise ValueError('Every public command must have exactly one workflow category')
    specialized=[]
    for row in SPECIALIZED:
        group=SPECIALIZED_CATEGORIES[row['kind']]
        command=('python '+row['script'] if row.get('script') else 'python -m '+row['module'])+' --help'
        specialized.append(dict(row,category=group,command=command,
            maintenance='compatibility_alias' if group=='historical' else 'maintained'))
    commands=[dict(command='redox '+name, category=group,
        module=COMMANDS[name][0], description=COMMANDS[name][2])
        for group,names in CATEGORIES.items() for name in names if name in COMMANDS]
    if category is not None:
        if category not in CATEGORIES:raise ValueError('Unknown workflow category')
        commands=[r for r in commands if r['category']==category]
        specialized=[r for r in specialized if r['category']==category]
    return dict(commands=commands,
        specialized=specialized, default='redox motifs plans; --execute is required for submission',
        qualification='A maintained entry point does not qualify every chemistry, GPU class or protocol')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--json',action='store_true')
    parser.add_argument('--category',choices=tuple(CATEGORIES),help='show only this maintained role; imports/submission are not performed')
    args=parser.parse_args(argv);data=catalogue(args.category)
    if args.json:print(json.dumps(data,indent=2));return
    for category in CATEGORIES:
        commands=[r for r in data['commands'] if r['category']==category]
        specialized=[r for r in data['specialized'] if r['category']==category]
        if not commands and not specialized:continue
        print(category.replace('_',' ').upper())
        for row in commands:print('  '+row['command'].ljust(28)+row['description'])
        for row in specialized:print('  '+row['command']+'\n    '+row['purpose'])
    print(data['default']);print(data['qualification'])


if __name__=='__main__':main()
