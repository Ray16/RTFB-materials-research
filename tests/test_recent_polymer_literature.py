import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).resolve().parents[1]
definition=importlib.util.spec_from_file_location('recent_polymer_citations',ROOT/'scripts/validation/check_recent_polymer_literature.py')
checker=importlib.util.module_from_spec(definition);definition.loader.exec_module(checker)


@pytest.fixture
def bibliography(tmp_path):
    source=ROOT/'validation/route_evidence_followup_2026-10-30'
    for name in ('recent_references.json','styrofoam_crossref.json','phthalimide_crossref.json'):
        shutil.copyfile(source/name,tmp_path/name)
    bib=tmp_path/'references.bib';shutil.copyfile(ROOT/'manuscript/novelty_references.bib',bib)
    return bib,tmp_path/'recent_references.json'


def test_both_recent_references_match_primary_records(bibliography):
    result=checker.audit(*bibliography)
    assert result['references']==2 and not result['route_or_property_promotion_inferred']


def test_semantic_comparison_accepts_line_wrapping_and_protected_capitalization(bibliography):
    bib,sources=bibliography
    text=bib.read_text().replace('Upcycling Commodity Styrofoam','Upcycling\n {Commodity Styrofoam}')
    text=text.replace(' and Sabrina',' and\n Sabrina').replace(' and Trano',' and\n Trano')
    bib.write_text(text)
    assert checker.audit(bib,sources)['references']==2


@pytest.mark.parametrize('mutation',['doi','author','pages','duplicate','metadata'])
def test_wrong_reference_identity_or_primary_record_is_refused(bibliography,mutation):
    bib,sources=bibliography;text=bib.read_text()
    if mutation=='doi':text=text.replace('cssc.202500092','cssc.000000000')
    if mutation=='author':text=text.replace('Oh Biying, Algin','Another, Author')
    if mutation=='pages':text=text.replace('9530--9541','9530--9999')
    if mutation=='duplicate':text+='\n'+text
    if mutation=='metadata':
        p=sources.parent/'styrofoam_crossref.json';p.write_text('{}')
    bib.write_text(text)
    with pytest.raises(ValueError):checker.audit(bib,sources)
