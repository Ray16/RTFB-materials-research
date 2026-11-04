"""A primary supplier graph supports exactly its positional N-H precursor."""
import json
import pytest
from redox.core.common import ROOT
from redox.workflows.molecular_native_namespace import load


def inputs():
    root=ROOT/'validation/phenothiazine_precursor_access_2026-11-02'
    panel=json.loads((ROOT/'results/candidates/pool/2026-11-02_direct_phenothiazine/inputs/ptz_direct_ring_2026-11-02.json').read_text())
    return panel,json.loads((root/'precursor_properties.json').read_text()),json.loads((root/'supplier_sid_cids.json').read_text())


def test_listing_matches_only_the_exact_2_cyano_NH_precursor():
    api=load(ROOT/'scripts/validation/check_ptz_precursor_access.py','tested_ptz_access')
    matches=api.matching_precursors(*inputs())
    assert [r['id'] for r in matches]==['proposal_ptz_176c4d2f397734ea']
    assert matches[0]['properties_computed'] is False


@pytest.mark.parametrize('bad',['supplier','CID','structure_key'])
def test_supplier_identity_or_graph_mismatch_is_rejected(bad):
    api=load(ROOT/'scripts/validation/check_ptz_precursor_access.py','tested_ptz_access')
    panel,properties,supplier=inputs()
    if bad=='supplier':supplier['InformationList']['Information'][0]['SID']=1
    if bad=='CID':supplier['InformationList']['Information'][0]['CID']=[374035]
    if bad=='structure_key':properties['PropertyTable']['Properties'][0]['InChIKey']='another positional isomer'
    with pytest.raises(ValueError):api.matching_precursors(panel,properties,supplier)
