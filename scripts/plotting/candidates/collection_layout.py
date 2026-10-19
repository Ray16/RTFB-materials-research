"""Stable, explicit destinations for the maintained candidate collection."""
from pathlib import PurePosixPath
import re


def artifact_path(name):
    name=str(name)
    if '/' in name:return name
    if name in ('index.html','README.md','.gitignore'):return name
    stem=PurePosixPath(name).stem
    suffix=PurePosixPath(name).suffix
    if name in ('manifest.json','captions.md','candidate_comparison_collection.tex'):
        return 'metadata/'+name
    if suffix in ('.json','.csv'):return 'data/'+name
    if stem in ('redox_capacity','synthesis_tradeoff','new_candidate_shifts','candidate_selection_overview'):
        return 'plots/'+name
    if stem in ('coverage','redox_ladders'):return 'supplementary/'+name
    if stem in ('candidate_structures','core_topology_comparison') or stem.startswith('redox_states'):
        return 'structures/'+name
    raise ValueError('Declare an output category before publishing: '+name)


def route_html(page,filenames):
    assets=set(filenames)
    return re.sub(r'\b(href|src)="([^"]+)"',
        lambda m:f'{m[1]}="{artifact_path(m[2])}"' if m[2] in assets else m[0],page)


def organize_stage(stage):
    files=[p for p in stage.iterdir() if p.is_file()]
    html=stage/'index.html'
    html.write_text(route_html(html.read_text(),[p.name for p in files]+['manifest.json']))
    for path in files:
        target=stage/artifact_path(path.name)
        if target!=path:
            target.parent.mkdir(parents=True,exist_ok=True);path.rename(target)
