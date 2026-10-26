#!/usr/bin/env python
"""Prepare a source-bound NDI external challenge without scoring a benchmark."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from xml.etree import ElementTree

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


SMILES = 'C=Cc1ccc(N2C(=O)c3ccc4c5c(ccc(c35)C2=O)C(=O)N(c2ccc(C=C)cc2)C4=O)cc1'
DOI = '10.3390/ijms26104807'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(xml_path, library_path, search_path):
    root = ElementTree.parse(xml_path).getroot()
    dois = [node.text for node in root.iter('article-id') if node.get('pub-id-type') == 'doi']
    if not dois or dois[0] != DOI:
        raise ValueError('Unexpected primary-source DOI')
    paragraphs = [' '.join(''.join(node.itertext()).split()) for node in root.iter('p')]
    chemistry = next((p for p in paragraphs if 'The targeted styrene–NDI monomer' in p), '')
    experiment = next((p for p in paragraphs if '−1040 mV' in p and '−1440 mV' in p), '')
    synthesis = next((p for p in paragraphs if 'C30H18N2O4' in p and '470.1267' in p), '')
    if not all(token in chemistry for token in ('4-vinylaniline', '1,4,5,8-naphthalenetetracarboxylic dianhydride')):
        raise ValueError('Expected source synthesis and starting materials')
    if not all(token in experiment for token in ('−880 mV', '−1350 mV', 'ferrocene/ferrocenium', 'acetonitrile', 'TBAClO4', 'quasi-reversible')):
        raise ValueError('Expected source electrode, solvent and peak descriptions')
    if not synthesis:
        raise ValueError('Expected source formula and mass characterization')
    mol = Chem.MolFromSmiles(SMILES)
    if mol is None or rdMolDescriptors.CalcMolFormula(mol) != 'C30H18N2O4':
        raise ValueError('Candidate graph does not match reported formula')
    checks = {'N_para_vinylphenyl': '[N](-c1ccc([CH]=[CH2])cc1)',
              'imide_N': '[N]([C]=[O])[C]=[O]', 'terminal_vinyl': '[CH]=[CH2]'}
    graph_counts = {name: len(mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)))
                    for name, pattern in checks.items()}
    if any(value != 2 for value in graph_counts.values()):
        raise ValueError('Graph lacks the two para-vinylphenyl imide arms')
    exact_mass = rdMolDescriptors.CalcExactMolWt(mol)
    if abs(exact_mass - 470.1267) > 0.0001:
        raise ValueError('Candidate graph does not match reported exact mass')
    key = Chem.MolToInchiKey(mol)
    library = list(csv.DictReader(Path(library_path).open(newline='')))
    library_matches = sorted({row['id'] for row in library if (m := Chem.MolFromSmiles(row['smiles'])) and Chem.MolToInchiKey(m) == key})
    search = list(csv.DictReader(Path(search_path).open(newline='')))
    search_matches = [row['id'] for row in search if row['inchikey'] == key]
    peaks = [dict(wave=1, cathodic_V=-1.040, anodic_V=-0.880),
             dict(wave=2, cathodic_V=-1.440, anodic_V=-1.350)]
    for row in peaks:
        row['peak_midpoint_V'] = (row['cathodic_V'] + row['anodic_V']) / 2
        row['peak_separation_V'] = row['anodic_V'] - row['cathodic_V']
    return dict(status='source_checked_NDI_external_challenge_not_scored', doi=DOI,
        source_url='https://doi.org/' + DOI, source_xml_sha256=digest(xml_path),
        molecule=dict(name='N,N-bis(4-vinylphenyl)-1,4,5,8-naphthalenediimide',
            smiles=Chem.MolToSmiles(mol), inchikey=key, formula=rdMolDescriptors.CalcMolFormula(mol),
            exact_mass_Da=exact_mass, graph_checks=graph_counts, source_figures=['Figure 1','Figure 5A'],
            identity_basis='source synthesis, formula/mass and visual check of Figure 1; proposed exact graph'),
        measurement=dict(solvent='acetonitrile', supporting_electrolyte='0.1 M tetrabutylammonium perchlorate',
            reference='external Fc/Fc+ calibration', scan_rate_monomer='not unambiguously stated in figure caption',
            temperature='not explicitly reported for monomer CV', behavior='quasi-reversible', peaks=peaks,
            midpoint_scope='Arithmetic peak midpoints are derived proxies, not source-reported formal potentials or validated equilibrium values.'),
        proposed_states=[dict(label='neutral',charge=0,multiplicity=1),dict(label='reduced_1',charge=-1,multiplicity=2),
                         dict(label='reduced_2',charge=-2,multiplicity=1)],
        state_scope='Parity-compatible hypotheses; alternative spin states, geometry and native thermochemistry unverified.',
        registry_comparison=dict(library_matches=library_matches, systematic_search_matches=search_matches,
            library_sha256=digest(library_path), systematic_search_sha256=digest(search_path)),
        role='Free-molecule experimental-method control only; not an exact graft or synthesis proposal.',
        strict_benchmark_admission=False, strict_benchmark_reason='Peak-based quasireversible midpoints and exact calculated states need review; fixed retrospective panel remains unchanged.',
        quantum_jobs_submitted=0, computed_potentials=None, experimental_accuracy_inferred=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--xml', type=Path, required=True)
    p.add_argument('--library', type=Path, required=True)
    p.add_argument('--search', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    result = prepare(args.xml, args.library, args.search)
    result['producer_sha256'] = digest(__file__)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
