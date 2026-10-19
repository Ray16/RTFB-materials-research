"""A withdrawn field cannot leak into the molecular comparison or its insight."""
import importlib.util
from pathlib import Path
import sys
import json
import pytest


def test_withdrawn_potential_is_blank_in_graph_plate_and_comparison(tmp_path):
    pytest.importorskip('cairosvg')
    folder=Path(__file__).resolve().parents[1]/'scripts/plotting/candidates';sys.path.insert(0,str(folder))
    from core_topology_plate import core_topology_plate
    smiles=['CCN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(Cc1ccc(C)cc1)C3=O','CCn1c(=O)c2c3ccccc3c3c(=O)n(Cc4ccc(C)cc4)c(=O)c3c2c1=O']
    rows=[dict(id=gid,exact_grafted_smiles=graph,E1_vs_Fc_V=e1,E2_vs_Fc_V=e2,capacity_nominal_mAh_g=134.5,potential_display_status=status) for gid,graph,e1,e2,status in zip(('design_ndi_ethyl','hyp_cata_ndi_ethyl'),smiles,(-.764,-.984),(-1.418,-1.730),('Pending / excluded','Admitted calculation'))]
    report=core_topology_plate(rows,tmp_path)
    assert report['same_formula'] and not report['comparison_published'] and report['delta']=={}
    assert report['structures'][0]['E1_vs_Fc_V'] is None
    svg=(tmp_path/'core_topology_comparison.svg').read_text()
    assert 'E1: pending' in svg and 'E1: -0.764' not in svg
    assert 'ΔE₁ =' not in svg
    assert json.loads((tmp_path/'core_topology_comparison.json').read_text())==report
