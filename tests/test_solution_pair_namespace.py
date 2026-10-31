from pathlib import Path
import shutil

import pytest

from redox.workflows.solution_pair_namespace import geometry_pair,command,validate_layout

ROOT=Path(__file__).resolve().parents[1]


def test_identical_real_cached_geometry_is_a_same_basin_pair():
    path=ROOT/'examples/cached_thermochemistry/data/product_neutral/opt.xyz'
    result=geometry_pair([path,path])
    assert result['heavy_RMSD_A']<1e-6
    assert result['maximum_bridge_torsion_difference_deg']<1e-6
    assert not result['global_conformer_convergence_proven']


def test_changed_atom_composition_cannot_be_used_as_a_shared_AO_pair(tmp_path):
    source=ROOT/'examples/cached_thermochemistry/data/product_neutral/opt.xyz'
    target=tmp_path/'different.xyz';lines=source.read_text().splitlines();fields=lines[2].split();fields[0]='N';lines[2]=' '.join(fields);target.write_text('\n'.join(lines)+'\n')
    with pytest.raises(ValueError):geometry_pair([source,target])


def test_two_matching_geometries_still_must_match_the_declared_candidate():
    path=ROOT/'examples/cached_thermochemistry/data/product_neutral/opt.xyz'
    with pytest.raises(ValueError,match='declared full graph'):
        geometry_pair([path,path],'CCO')


def layout(tmp_path):
    return dict(prior_root='/prior',archive='/archive',case=str(tmp_path/'case'),sources=str(tmp_path/'sources'))


def test_source_tree_is_read_only_and_only_declared_case_is_writable(tmp_path):
    spec=layout(tmp_path);argv=command(spec)
    assert argv[argv.index('--ro-bind',argv.index('--ro-bind')+1)+1:][:2]==['/archive','/prior']
    assert argv[argv.index('--bind')+1:argv.index('--bind')+3]==[spec['case'],spec['case']]
    assert '--tmpfs' in argv


@pytest.mark.parametrize('change',[dict(case='/prior/results'),dict(case='/archive/results'),
                                  dict(sources='/archive/code'),dict(case='/'),dict(case='relative'),
                                  dict(case='/other/../case'),dict(sources='/another-root/sources')])
def test_unsafe_or_ambiguous_writable_mounts_are_refused(tmp_path,change):
    spec=layout(tmp_path);spec.update(change)
    with pytest.raises(ValueError):validate_layout(spec)
