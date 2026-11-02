"""Read-only catalogue of maintained commands, controls and recovery entry points."""
import argparse
import json

CATEGORIES = {
 'production': ['motifs','launch','run-gpu'],
 'design': ['build','build-validation','expand','search','source-audit'],
 'analysis_and_validation': ['status','energy-audit','review','thermal-postprocess','uma-check','seed-reference','resources'],
 'calculation_stages': ['uma','dft','sp','native-thermal'],
 'discovery': ['workflows'],
}
SPECIALIZED = [
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


def catalogue():
    from redox.cli import COMMANDS
    listed = [name for group in CATEGORIES.values() for name in group if name in COMMANDS]
    if len(set(listed)) != len(listed) or set(listed) != set(COMMANDS):
        raise ValueError('Every public command must have exactly one workflow category')
    return dict(commands=[dict(command='redox '+name, category=category,
        module=COMMANDS[name][0], description=COMMANDS[name][2])
        for category,names in CATEGORIES.items() for name in names if name in COMMANDS],
        specialized=SPECIALIZED, default='redox motifs plans; --execute is required for submission',
        qualification='A maintained entry point does not qualify every chemistry, GPU class or protocol')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--json',action='store_true');args=parser.parse_args(argv)
    data=catalogue()
    if args.json:print(json.dumps(data,indent=2));return
    for category in CATEGORIES:
        print(category.replace('_',' ').upper())
        for row in data['commands']:
            if row['category']==category:print('  '+row['command'].ljust(28)+row['description'])
    print('SPECIALIZED CAMPAIGNS / CONTROLS / RECOVERY')
    for row in data['specialized']:print('  '+row['kind']+': python -m '+row['module']+' --help\n    '+row['purpose'])
    print(data['default']);print(data['qualification'])


if __name__=='__main__':main()
