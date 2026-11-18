"""Registered grafted motifs. Provenance is historical; current property QC comes from addressed calculation records.

Exactly one terminal [*:1] attachment and explicit (label, charge, multiplicity, n_e) ladder per motif.
Source parent measurements are not properties of these grafts.
"""
SCAFFOLD = 'Cc1ccc(C[*:1])cc1'
GROUPS = [{'id': 'exp_viologen_methoxyethyl',
  'name': "N-(4-methylbenzyl)-N'-(2-methoxyethyl)-4,4'-bipyridinium",
  'family': 'pyridine-multi-e',
  'frag': '[*:1][n+]1ccc(-c2cc[n+](CCOC)cc2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'provenance': {'type': 'uncomputed grafted proposal',
                 'source': 'https://patents.google.com/patent/US8773745B2/en',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'exact ether-cap precursor and benzylic-bromide polymer graft '
                                 'published; chloride transfer proposed',
                 'direct_merrifield_route_verified': False,
                 'experimental_parent_reference': '10.1016/j.jechem.2018.02.014',
                 'registration_stage': 'grafting and UMA warm-start only',
                 'registered': '2026-10-11'},
  'chemistry_notes': 'US8773745B2 publishes the exact monoquaternized ether-cap precursor and its '
                     'attachment to benzylic-bromide polystyrene. Borchers2022 demonstrates '
                     'chloromethylstyrene attachment with a different cap. Ether-cap transfer to '
                     'Merrifield chloride remains to establish.',
  'review_priority': 'pilot; physical metric gate pending'},
 {'id': 'exp_phthalimide_control',
  'name': 'N-(4-methylbenzyl)phthalimide',
  'family': 'mono-imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccccc2C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1]],
  'provenance': {'type': 'uncomputed grafted proposal',
                 'source': 'https://pubs.rsc.org/en/content/articlehtml/2019/ra/c9ra01131e',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'direct phthalimide substitution onto chloromethyl-polystyrene '
                                 'published; exact truncated proxy properties uncomputed',
                 'direct_merrifield_route_verified': True,
                 'experimental_parent_reference': '10.1021/acsaem.1c01490',
                 'registration_stage': 'grafting and UMA warm-start only',
                 'registered': '2026-10-11'},
  'chemistry_notes': 'RSC2019 C9RA01131E reports potassium phthalimide substitution onto '
                     'poly(4-chloromethylstyrene). This flag means a published direct support-reaction '
                     'precedent, not experimental execution here or validated electrochemical '
                     'performance. Daub provides monoimide one-electron evidence.',
  'review_priority': 'pilot; physical metric gate pending'},
 {'id': 'exp_ndi_hexyl_pentyloxy',
  'name': "N-[5-(4-methylbenzyloxy)pentyl]-N'-n-hexyl naphthalene diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]OCCCCCN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCCCCC)C3=O',
  'states': [('neu', 0, 1, 0), ('red1', -1, 2, -1), ('red2', -2, 1, -2)],
  'provenance': {'type': 'uncomputed grafted proposal',
                 'source': 'https://pubs.rsc.org/en/content/articlehtml/2021/qm/d0qm00685h',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'exact isolated hydroxyalkyl precursor published; terminal-alcohol '
                                 'chloride graft proposed',
                 'direct_merrifield_route_verified': False,
                 'experimental_parent_reference': '10.1039/D0QM00685H',
                 'registration_stage': 'grafted inputs prepared; solution properties pending',
                 'registered': '2026-10-11'},
  'chemistry_notes': 'This removes the unestablished mono-N-H diimide requirement. The paper '
                     'esterifies the OH with norbornenecarboxylic acid; it does not report this ether '
                     'graft. Its redox data are in DCM and are not MeCN property anchors.',
  'review_priority': 'balanced first calculation batch; physical protocol gate pending'},
 {'id': 'exp_aq_1oxy',
  'name': '1-(4-methylbenzyloxy)anthracene-9,10-dione',
  'family': 'quinone (n-type)',
  'frag': '[*:1]Oc1cccc2c1C(=O)c1ccccc1C2=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'uncomputed grafted proposal',
                 'source': 'https://patents.google.com/patent/WO2004084891A1/en',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'hydroxy-AQ benzyl ether reaction-class precedent; exact 1-O resin '
                                 'route proposed',
                 'direct_merrifield_route_verified': False,
                 'experimental_parent_reference': '10.1021/acs.jpcc.2c03129',
                 'registration_stage': 'grafted inputs prepared; solution properties pending',
                 'registered': '2026-10-11'},
  'chemistry_notes': 'The cited reaction analog has a nitro substituent absent from this target. 1-O '
                     'is peri to a carbonyl; this changes intramolecular H-bonding relative to parent. '
                     'Parent benchmark half-peak potentials are not exact thermodynamic E values.',
  'review_priority': 'balanced first calculation batch; physical protocol gate pending'},
 {'id': 'design_pmdi_ethyl',
  'name': "N-(4-methylbenzyl)-N'-ethyl pyromellitic diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n(CC)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'series',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'design_pmdi_methoxyethyl',
  'name': "N-(4-methylbenzyl)-N'-(2-methoxyethyl) pyromellitic diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n(CCOC)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'next chemistry smoke',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'design_pmdi_methoxyethoxyethyl',
  'name': "N-(4-methylbenzyl)-N'-[2-(2-methoxyethoxy)ethyl] pyromellitic diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n(CCOCCOC)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'series',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'design_pmdi_isopropyl',
  'name': "N-(4-methylbenzyl)-N'-isopropyl pyromellitic diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n(C(C)C)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'series',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'design_ndi_methyl',
  'name': "N-(4-methylbenzyl)-N'-methyl naphthalene diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(C)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['ndi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'next chemistry smoke',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'design_ndi_ethyl',
  'name': "N-(4-methylbenzyl)-N'-ethyl naphthalene diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CC)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['ndi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'series',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'design_ndi_methoxyethyl',
  'name': "N-(4-methylbenzyl)-N'-(2-methoxyethyl) naphthalene diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCOC)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['ndi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'next chemistry smoke',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'design_ndi_methoxyethoxyethyl',
  'name': "N-(4-methylbenzyl)-N'-[2-(2-methoxyethoxy)ethyl] naphthalene diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCOCCOC)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['ndi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'series',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'design_ndi_isopropyl',
  'name': "N-(4-methylbenzyl)-N'-isopropyl naphthalene diimide",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(C(C)C)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['ndi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior nine-imide proposal review',
                 'source': 'https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'review_priority': 'series',
  'chemistry_notes': 'Single imide N-H precursor must be selectively prepared and isolated; '
                     'indiscriminate mixed-amine dianhydride condensation can give '
                     'symmetric/asymmetric mixtures. Both imide units/four carbonyls retained. Other N '
                     'capped; no residual imide NH. Parent paper supports scaffold, not measurements '
                     'of this benzyl structure.'},
 {'id': 'd3_05ujsa',
  'name': 'D3TaLES 05UJSA N-benzyl (amine N-H) product (canonical structure supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CCCCN(C1=C(Cl)C(=O)c2ccccc2C1=O)[*:1]',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '05UJSA',
                 'parent_smiles': 'CCCCNC1=C(Cl)C(=O)c2ccccc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.4833',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Ring halogen: competing nucleophilic substitution must be checked. Conjugated '
                     'aminoquinone N is a weak nucleophile; modeled graft is not a verified direct '
                     'Merrifield substitution route. Resin aminomethyl addition to quinone is an '
                     'alternate architecture. D3TaLES lambda is parent computed evidence only; no '
                     'second-step grafted lambda or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_05lcla',
  'name': 'D3TaLES 05LCLA benzyl ester (carboxylate) product (canonical structure supplied)',
  'family': 'quinone (n-type)',
  'frag': 'O=C(CCSC1=CC(=O)C(=O)c2ccccc21)O[*:1]',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '05LCLA',
                 'parent_smiles': 'O=C(O)CCSC1=CC(=O)C(=O)c2ccccc21',
                 'd3_parent_electron_reorg_eV': '0.6269',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Ester linkage cleavage and reactive o-quinone chemistry require stability '
                     'review. D3TaLES lambda is parent computed evidence only; no second-step grafted '
                     'lambda or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_05mnmt',
  'name': 'D3TaLES 05MNMT N-benzyl (amine N-H) product (canonical structure supplied)',
  'family': 'quinone (n-type)',
  'frag': 'O=C1C(Br)=C(N[*:1])C(=O)c2ccccc21',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '05MNMT',
                 'parent_smiles': 'NC1=C(Br)C(=O)c2ccccc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.507',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Ring halogen: competing nucleophilic substitution must be checked. Conjugated '
                     'aminoquinone N is a weak nucleophile; modeled graft is not a verified direct '
                     'Merrifield substitution route. Resin aminomethyl addition to quinone is an '
                     'alternate architecture. Residual OH/NH donor remains; proton transfer and '
                     'competitive attachment require review. D3TaLES lambda is parent computed '
                     'evidence only; no second-step grafted lambda or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_05txcn',
  'name': 'D3TaLES 05TXCN N-benzyl (amine N-H) product (canonical structure supplied)',
  'family': 'quinone (n-type)',
  'frag': 'O=C1C(Cl)=C(N[*:1])C(=O)c2ccccc21',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '05TXCN',
                 'parent_smiles': 'NC1=C(Cl)C(=O)c2ccccc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.5506',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Ring halogen: competing nucleophilic substitution must be checked. Conjugated '
                     'aminoquinone N is a weak nucleophile; modeled graft is not a verified direct '
                     'Merrifield substitution route. Resin aminomethyl addition to quinone is an '
                     'alternate architecture. Residual OH/NH donor remains; proton transfer and '
                     'competitive attachment require review. D3TaLES lambda is parent computed '
                     'evidence only; no second-step grafted lambda or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_05czhi',
  'name': 'D3TaLES 05CZHI O-benzyl ether (phenol/enol, Williamson) product (canonical structure '
          'supplied)',
  'family': 'quinone (n-type)',
  'frag': 'O=C1C(Cl)=C(Cl)C(=O)c2c(O[*:1])ccc(O)c21',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '05CZHI',
                 'parent_smiles': 'O=C1C(Cl)=C(Cl)C(=O)c2c(O)ccc(O)c21',
                 'd3_parent_electron_reorg_eV': '0.2449',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Ring halogen: competing nucleophilic substitution must be checked. Peripheral '
                     'phenol/alcohol O-benzylation leaves two quinone carbonyls intact. Residual OH/NH '
                     'donor remains; proton transfer and competitive attachment require review. '
                     'D3TaLES lambda is parent computed evidence only; no second-step grafted lambda '
                     'or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_80xwvs',
  'name': 'D3TaLES 80XWVS N-benzyl (amine N-H) product (canonical structure supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CN(C1=CC(=O)c2ccncc2C1=O)[*:1]',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '80XWVS',
                 'parent_smiles': 'CNC1=CC(=O)c2ccncc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.6123',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Conjugated aminoquinone N is a weak nucleophile; modeled graft is not a verified '
                     'direct Merrifield substitution route. Resin aminomethyl addition to quinone is '
                     'an alternate architecture. Aza ring N competes through Menshutkin alkylation and '
                     'can change core charge; attachment selectivity required. D3TaLES lambda is '
                     'parent computed evidence only; no second-step grafted lambda or potential is '
                     'known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_80mdpm',
  'name': 'D3TaLES 80MDPM N-benzyl (amine N-H) product (canonical structure supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CN(C1=C(Cl)C(=O)c2ccncc2C1=O)[*:1]',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '80MDPM',
                 'parent_smiles': 'CNC1=C(Cl)C(=O)c2ccncc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.4958',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Ring halogen: competing nucleophilic substitution must be checked. Conjugated '
                     'aminoquinone N is a weak nucleophile; modeled graft is not a verified direct '
                     'Merrifield substitution route. Resin aminomethyl addition to quinone is an '
                     'alternate architecture. Aza ring N competes through Menshutkin alkylation and '
                     'can change core charge; attachment selectivity required. D3TaLES lambda is '
                     'parent computed evidence only; no second-step grafted lambda or potential is '
                     'known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_80necn',
  'name': 'D3TaLES 80NECN O-benzyl ether (phenol/enol, Williamson) product (canonical structure '
          'supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CN(C)C1=CC(=O)c2cc(O[*:1])ccc2C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'additional D3TaLES pool candidate',
                 'source_id': '80NECN',
                 'parent_smiles': 'CN(C)C1=CC(=O)c2cc(O)ccc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.6776',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Peripheral phenol/alcohol O-benzylation leaves two quinone carbonyls intact. '
                     'Conjugated tertiary amine may compete by quaternization; single phenolic handle '
                     'is not proof of selectivity. D3TaLES lambda is parent computed evidence only; no '
                     'second-step grafted lambda or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_80nopt',
  'name': 'D3TaLES 80NOPT N-benzyl (amine N-H) product (canonical structure supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CCCCN(C1=CC(=O)c2ccccc2C1=O)[*:1]',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '80NOPT',
                 'parent_smiles': 'CCCCNC1=CC(=O)c2ccccc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.5311',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Conjugated aminoquinone N is a weak nucleophile; modeled graft is not a verified '
                     'direct Merrifield substitution route. Resin aminomethyl addition to quinone is '
                     'an alternate architecture. D3TaLES lambda is parent computed evidence only; no '
                     'second-step grafted lambda or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_80levd',
  'name': 'D3TaLES 80LEVD O-benzyl ether (phenol/enol, Williamson) product (canonical structure '
          'supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CC(C)=CCC1=C(O[*:1])C(=O)c2ccccc2C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '80LEVD',
                 'parent_smiles': 'CC(C)=CCC1=C(O)C(=O)c2ccccc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.4171',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Enolic OH: O/C alkylation selectivity required; both carbonyls must survive. '
                     'D3TaLES lambda is parent computed evidence only; no second-step grafted lambda '
                     'or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_80prff',
  'name': 'D3TaLES 80PRFF O-benzyl ether (phenol/enol, Williamson) product (canonical structure '
          'supplied)',
  'family': 'quinone (n-type)',
  'frag': 'O=C1C=CC(=O)c2c(O[*:1])cccc21',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'additional D3TaLES pool candidate',
                 'source_id': '80PRFF',
                 'parent_smiles': 'O=C1C=CC(=O)c2c(O)cccc21',
                 'd3_parent_electron_reorg_eV': '0.4552',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Peripheral phenol/alcohol O-benzylation leaves two quinone carbonyls intact. '
                     'Peri phenolic O graft removes OH-to-carbonyl H-bond; recompute grafted '
                     'properties. D3TaLES lambda is parent computed evidence only; no second-step '
                     'grafted lambda or potential is known.',
  'review_priority': 'next chemistry smoke'},
 {'id': 'd3_80gkvh',
  'name': 'D3TaLES 80GKVH O-benzyl ether (phenol/enol, Williamson) product (canonical structure '
          'supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CC(C)=CC1=C(O[*:1])C(=O)c2ccccc2C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'prior ten-front review',
                 'source_id': '80GKVH',
                 'parent_smiles': 'CC(C)=CC1=C(O)C(=O)c2ccccc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.4671',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Enolic OH: O/C alkylation selectivity required; both carbonyls must survive. '
                     'D3TaLES lambda is parent computed evidence only; no second-step grafted lambda '
                     'or potential is known.',
  'review_priority': 'defer route review'},
 {'id': 'd3_80gxwg',
  'name': 'D3TaLES 80GXWG O-benzyl ether (phenol/enol, Williamson) product (canonical structure '
          'supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CC1=C(CO[*:1])C(=O)c2ccc(C)cc2C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'additional D3TaLES pool candidate',
                 'source_id': '80GXWG',
                 'parent_smiles': 'CC1=C(CO)C(=O)c2ccc(C)cc2C1=O',
                 'd3_parent_electron_reorg_eV': '1.558',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Peripheral phenol/alcohol O-benzylation leaves two quinone carbonyls intact. '
                     'D3TaLES lambda is parent computed evidence only; no second-step grafted lambda '
                     'or potential is known.',
  'review_priority': 'series'},
 {'id': 'd3_80eeoz',
  'name': 'D3TaLES 80EEOZ O-benzyl ether (phenol/enol, Williamson) product (canonical structure '
          'supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CC1=C(CCCO[*:1])C(=O)c2ccc(C)cc2C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'additional D3TaLES pool candidate',
                 'source_id': '80EEOZ',
                 'parent_smiles': 'CC1=C(CCCO)C(=O)c2ccc(C)cc2C1=O',
                 'd3_parent_electron_reorg_eV': '1.411',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Peripheral phenol/alcohol O-benzylation leaves two quinone carbonyls intact. '
                     'D3TaLES lambda is parent computed evidence only; no second-step grafted lambda '
                     'or potential is known.',
  'review_priority': 'next chemistry smoke'},
 {'id': 'd3_80tbwr',
  'name': 'D3TaLES 80TBWR O-benzyl ether (phenol/enol, Williamson) product (canonical structure '
          'supplied)',
  'family': 'quinone (n-type)',
  'frag': 'CC1=C(CCO[*:1])C(=O)c2occc2C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['d3tales', 'quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'additional D3TaLES pool candidate',
                 'source_id': '80TBWR',
                 'parent_smiles': 'CC1=C(CCO)C(=O)c2occc2C1=O',
                 'd3_parent_electron_reorg_eV': '0.6437',
                 'source': 'https://doi.org/10.18126/v5sj-6q93; '
                           'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Peripheral phenol/alcohol O-benzylation leaves two quinone carbonyls intact. '
                     'D3TaLES lambda is parent computed evidence only; no second-step grafted lambda '
                     'or potential is known.',
  'review_priority': 'series'},
 {'id': 'exp_pmdi_n_butyl',
  'name': "N-(4-methylbenzyl)-N'-butyl PMDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n(CCCC)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Exactly one imide N-H precursor; isolate asymmetric mono-capped precursor before '
                     'N-benzylation. Both imide units/four carbonyls retained. No residual donor '
                     'proton. Longer chains add mass and conformers; no assumed solubility advantage '
                     'in polymer phase.',
  'review_priority': 'series'},
 {'id': 'exp_pmdi_trioxyethyl',
  'name': "N-(4-methylbenzyl)-N'-2-[2-(2-methoxyethoxy)ethoxy]ethyl PMDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n(CCOCCOCCOC)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Exactly one imide N-H precursor; isolate asymmetric mono-capped precursor before '
                     'N-benzylation. Both imide units/four carbonyls retained. No residual donor '
                     'proton. Longer chains add mass and conformers; no assumed solubility advantage '
                     'in polymer phase.',
  'review_priority': 'series'},
 {'id': 'exp_ndi_n_butyl',
  'name': "N-(4-methylbenzyl)-N'-butyl NDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCCC)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['ndi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Exactly one imide N-H precursor; isolate asymmetric mono-capped precursor before '
                     'N-benzylation. Both imide units/four carbonyls retained. No residual donor '
                     'proton. Longer chains add mass and conformers; no assumed solubility advantage '
                     'in polymer phase.',
  'review_priority': 'series'},
 {'id': 'exp_ndi_trioxyethyl',
  'name': "N-(4-methylbenzyl)-N'-2-[2-(2-methoxyethoxy)ethoxy]ethyl NDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCOCCOCCOC)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['ndi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/cssc.202501011',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Exactly one imide N-H precursor; isolate asymmetric mono-capped precursor before '
                     'N-benzylation. Both imide units/four carbonyls retained. No residual donor '
                     'proton. Longer chains add mass and conformers; no assumed solubility advantage '
                     'in polymer phase.',
  'review_priority': 'series'},
 {'id': 'exp_pmdi_pent2_r',
  'name': "N-(4-methylbenzyl)-N'-(R)-pentan-2-yl PMDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n([C@H](C)CCC)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Explicit pentan-2-yl stereoisomer; CIP identity exported from RDKit. Source '
                     'parent uses branched pentyl substitution but proposed asymmetric tethered '
                     'enantiomer is unmeasured. Isolate mono-N-H precursor; four imide carbonyls '
                     'retained.',
  'review_priority': 'series'},
 {'id': 'exp_pmdi_pent2_s',
  'name': "N-(4-methylbenzyl)-N'-(S)-pentan-2-yl PMDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n([C@@H](C)CCC)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Explicit pentan-2-yl stereoisomer; CIP identity exported from RDKit. Source '
                     'parent uses branched pentyl substitution but proposed asymmetric tethered '
                     'enantiomer is unmeasured. Isolate mono-N-H precursor; four imide carbonyls '
                     'retained.',
  'review_priority': 'series'},
 {'id': 'exp_aq_2oxy_6methoxy',
  'name': '2-(4-methylbenzyloxy)-6-methoxyanthracene-9,10-dione',
  'family': 'quinone (n-type)',
  'frag': '[*:1]Oc1ccc2c(c1)C(=O)c1ccc(OC)cc1C2=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Single non-peri phenolic handle and capped second oxygen. Preserve quinone '
                     'carbonyls and avoid residual OH/PCET. Rational ether derivative; exact parent is '
                     'not claimed as a measured paper compound.',
  'review_priority': 'next chemistry smoke'},
 {'id': 'exp_bq_oxy_methoxy',
  'name': '2-(4-methylbenzyloxy)-5-methoxy-1,4-benzoquinone',
  'family': 'quinone (n-type)',
  'frag': '[*:1]OC1=CC(=O)C(OC)=CC1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Enolic hydroxyquinone precursor needs O/C alkylation selectivity check. Two '
                     'carbonyls retained and no residual OH. Unshielded ring positions may undergo '
                     'nucleophilic addition.',
  'review_priority': 'defer route review'},
 {'id': 'exp_bq_oxy_trimethyl',
  'name': '2-(4-methylbenzyloxy)-3,5,6-trimethyl-1,4-benzoquinone',
  'family': 'quinone (n-type)',
  'frag': '[*:1]OC1=C(C)C(=O)C(C)=C(C)C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'source_keys': ['quinone'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.acs.org/doi/10.1021/acs.jpcc.0c03632',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Rational methyl shielding series, not a verified measured compound. Enolic OH '
                     'precursor requires O/C selectivity review; methyl substitution does not prove '
                     'cycling stability.',
  'review_priority': 'defer route review'},
 {'id': 'exp_viologen_methoxyethoxyethyl',
  'name': "N-(4-methylbenzyl)-N'-[2-(2-methoxyethoxy)ethyl]-4,4'-bipyridinium",
  'family': 'pyridine-multi-e',
  'frag': '[*:1][n+]1ccc(-c2cc[n+](CCOCCOC)cc2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'provenance': {'source_keys': ['viologen', 'glycol_viologen'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://www.sciencedirect.com/science/article/pii/S2095495617311403; '
                           'https://pubs.rsc.org/en/content/articlelanding/2022/ta/d2ta07876g',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'One remaining pyridine N after mono-quaternization; +2/+1/0 ladder. Ether chain '
                     'is a deliberate polarity/mass control, not a solubility prediction.',
  'review_priority': 'series'},
 {'id': 'exp_naphthalimide_control',
  'name': 'N-(4-methylbenzyl)-1,8-naphthalimide',
  'family': 'mono-imide (n-type)',
  'frag': '[*:1]N1C(=O)c2cccc3cccc(c23)C1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1]],
  'provenance': {'source_keys': ['pmdi'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'Rational mono-imide core-size control, not a structure measured in cited PMDI '
                     'paper. One-electron proposed ladder; reversibility for this enlarged mono-imide '
                     'remains to test.',
  'review_priority': 'one-electron control'},
 {'id': 'exp_ptz_37dimethoxy',
  'name': 'N-(4-methylbenzyl)-3,7-dimethoxyphenothiazine',
  'family': 'phenothiazine (p-type)',
  'frag': '[*:1]N1c2ccc(OC)cc2Sc2cc(OC)ccc21',
  'states': [['neu', 0, 1, 0], ['ox1', 1, 2, 1]],
  'provenance': {'source_keys': ['phenothiazine', 'glycol_phenothiazine'],
                 'type': 'source-inspired grafted design',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'origin': 'new expansion',
                 'source': 'https://pubs.rsc.org/en/content/articlelanding/2017/ta/c7ta05883g; '
                           'https://pubs.acs.org/doi/10.1021/acs.chemmater.8b04770',
                 'registered': '2026-10-11',
                 'registration_stage': 'original39 exact graph/stereochemistry; computations pending',
                 'exact_grafted_experiment': False},
  'chemistry_notes': 'N-H alkylation preserves the two methoxy donors para to N that support the '
                     'parent radical cation. Reclassified 2026-11-16 to a single reversible 1e '
                     'oxidation (0/+1): the proposed second oxidation (+1/+2) was not reproduced by '
                     'the pipeline (benchmark error +0.55 V vs experiment on this dimethoxy core), '
                     'and phenothiazine is a validated 1e motif. Retained as a 1e catholyte '
                     'reference, not a two-electron candidate.',
  'review_priority': 'separate catholyte branch (1e)'}]


# Matched-comparison hypotheses registered 2026-10-18; exact-graft properties remain pending.
GROUPS += [{'id': 'hyp_pmdi_cyclopropyl',
  'name': "N-(4-methylbenzyl)-N'-cyclopropyl PMDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n(C5CC5)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'comparison_baseline': 'design_pmdi_isopropyl',
  'hypothesis': 'Test cap topology and conformer sensitivity against isopropyl; cyclization also '
                'removes two H atoms, so this is not a pure rigidity perturbation.',
  'provenance': {'type': 'matched-comparison computational hypothesis',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'evidence_scope': 'Family rationale for imide multiredox screening; no '
                                   'measurements or synthesis proof for these grafts.',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'direct_merrifield_route_verified': False,
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Test cap topology and conformer sensitivity against isopropyl; cyclization '
                     'also removes two H atoms, so this is not a pure rigidity perturbation. Exact '
                     'graft/precursor synthesis is unverified; parent measurements are not graft '
                     'properties. Imide mono-N-H precursor isolation or selective phenothiazine '
                     'N-alkylation requires separate route review.'},
 {'id': 'hyp_pmdi_2fluoroethyl',
  'name': "N-(4-methylbenzyl)-N'-2fluoroethyl PMDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]n1c(=O)c2cc3c(=O)n(CCF)c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'comparison_baseline': 'design_pmdi_ethyl',
  'hypothesis': 'Test replacing one terminal ethyl H by F against the ethyl cap; evaluate E1/E2, '
                'wave spacing and mass cost without assuming improved solubility.',
  'provenance': {'type': 'matched-comparison computational hypothesis',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'evidence_scope': 'Family rationale for imide multiredox screening; no '
                                   'measurements or synthesis proof for these grafts.',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'direct_merrifield_route_verified': False,
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Test replacing one terminal ethyl H by F against the ethyl cap; evaluate '
                     'E1/E2, wave spacing and mass cost without assuming improved solubility. '
                     'Exact graft/precursor synthesis is unverified; parent measurements are not '
                     'graft properties. Imide mono-N-H precursor isolation or selective '
                     'phenothiazine N-alkylation requires separate route review.'},
 {'id': 'hyp_ndi_cyclopropyl',
  'name': "N-(4-methylbenzyl)-N'-cyclopropyl NDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(C5CC5)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'comparison_baseline': 'design_ndi_isopropyl',
  'hypothesis': 'Test cap topology and conformer sensitivity against isopropyl; cyclization also '
                'removes two H atoms, so this is not a pure rigidity perturbation.',
  'provenance': {'type': 'matched-comparison computational hypothesis',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'evidence_scope': 'Family rationale for imide multiredox screening; no '
                                   'measurements or synthesis proof for these grafts.',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'direct_merrifield_route_verified': False,
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Test cap topology and conformer sensitivity against isopropyl; cyclization '
                     'also removes two H atoms, so this is not a pure rigidity perturbation. Exact '
                     'graft/precursor synthesis is unverified; parent measurements are not graft '
                     'properties. Imide mono-N-H precursor isolation or selective phenothiazine '
                     'N-alkylation requires separate route review.'},
 {'id': 'hyp_ndi_2fluoroethyl',
  'name': "N-(4-methylbenzyl)-N'-2fluoroethyl NDI",
  'family': 'imide (n-type)',
  'frag': '[*:1]N1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCF)C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'comparison_baseline': 'design_ndi_ethyl',
  'hypothesis': 'Test replacing one terminal ethyl H by F against the ethyl cap; evaluate E1/E2, '
                'wave spacing and mass cost without assuming improved solubility.',
  'provenance': {'type': 'matched-comparison computational hypothesis',
                 'source': 'https://pubs.acs.org/doi/10.1021/acsaem.1c01490',
                 'evidence_scope': 'Family rationale for imide multiredox screening; no '
                                   'measurements or synthesis proof for these grafts.',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'direct_merrifield_route_verified': False,
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Test replacing one terminal ethyl H by F against the ethyl cap; evaluate '
                     'E1/E2, wave spacing and mass cost without assuming improved solubility. '
                     'Exact graft/precursor synthesis is unverified; parent measurements are not '
                     'graft properties. Imide mono-N-H precursor isolation or selective '
                     'phenothiazine N-alkylation requires separate route review.'}]
# Removed 2026-11-16: hyp_ptz_37dimethyl and hyp_ptz_37difluoro. Both were speculative
# two-electron (0/+1/+2) phenothiazine proposals; phenothiazine is a validated 1e motif and
# the pipeline did not reproduce a second oxidation (benchmark +0.55 V error on the dimethoxy
# core). No production QM was spent on them. See exp_ptz_37dimethoxy, kept as a 1e catholyte.


# Matched-comparison hypotheses registered 2026-10-18; exact-graft properties remain pending.
GROUPS += [{'id': 'search_viologen_991d299c38be',
  'name': 'viologen / direct_N / npropyl / H',
  'family': 'pyridine-multi-e',
  'frag': 'CCC[n+]1ccc(-c2cc[n+]([*:1])cc2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'ethylviologen',
  'search_axes': {'core': 'viologen',
                  'anchor': 'direct_N',
                  'cap': 'npropyl',
                  'ring': 'H',
                  'template_atom': None},
  'comparison_contrast': 'cap_length',
  'hypothesis': 'Measure cap_length on viologen against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_viologen_1c2d9d5208a4',
  'name': 'viologen / direct_N / ethyl / F',
  'family': 'pyridine-multi-e',
  'frag': 'CC[n+]1ccc(-c2cc[n+]([*:1])c(F)c2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'ethylviologen',
  'search_axes': {'core': 'viologen',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'F',
                  'template_atom': 2},
  'comparison_contrast': 'ring_F',
  'hypothesis': 'Measure ring_F on viologen against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_viologen_e46c2caea23e',
  'name': 'viologen / direct_N / ethyl / CN',
  'family': 'pyridine-multi-e',
  'frag': 'CC[n+]1ccc(-c2cc[n+]([*:1])c(C#N)c2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'ethylviologen',
  'search_axes': {'core': 'viologen',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'CN',
                  'template_atom': 2},
  'comparison_contrast': 'ring_CN',
  'hypothesis': 'Measure ring_CN on viologen against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_viologen_6ddf8e94b88f',
  'name': 'viologen / ether_C2 / ethyl / H',
  'family': 'pyridine-multi-e',
  'frag': 'CC[n+]1ccc(-c2cc[n+](CCO[*:1])cc2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'ethylviologen',
  'search_axes': {'core': 'viologen',
                  'anchor': 'ether_C2',
                  'cap': 'ethyl',
                  'ring': 'H',
                  'template_atom': None},
  'comparison_contrast': 'linker_C2',
  'hypothesis': 'Measure linker_C2 on viologen against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_pmdi_0fe159524740',
  'name': 'pmdi / direct_N / npropyl / H',
  'family': 'imide (n-type)',
  'frag': 'CCCn1c(=O)c2cc3c(=O)n([*:1])c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_pmdi_ethyl',
  'search_axes': {'core': 'pmdi',
                  'anchor': 'direct_N',
                  'cap': 'npropyl',
                  'ring': 'H',
                  'template_atom': None},
  'comparison_contrast': 'cap_length',
  'hypothesis': 'Measure cap_length on pmdi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_pmdi_e5caac5870e1',
  'name': 'pmdi / direct_N / ethyl / F',
  'family': 'imide (n-type)',
  'frag': 'CCn1c(=O)c2cc3c(=O)n([*:1])c(=O)c3c(F)c2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_pmdi_ethyl',
  'search_axes': {'core': 'pmdi',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'F',
                  'template_atom': 5},
  'comparison_contrast': 'ring_F',
  'hypothesis': 'Measure ring_F on pmdi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_pmdi_1864c2c74267',
  'name': 'pmdi / direct_N / ethyl / CN',
  'family': 'imide (n-type)',
  'frag': 'CCn1c(=O)c2cc3c(=O)n([*:1])c(=O)c3c(C#N)c2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_pmdi_ethyl',
  'search_axes': {'core': 'pmdi',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'CN',
                  'template_atom': 5},
  'comparison_contrast': 'ring_CN',
  'hypothesis': 'Measure ring_CN on pmdi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_pmdi_87840c94153a',
  'name': 'pmdi / ether_C2 / ethyl / H',
  'family': 'imide (n-type)',
  'frag': 'CCn1c(=O)c2cc3c(=O)n(CCO[*:1])c(=O)c3cc2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_pmdi_ethyl',
  'search_axes': {'core': 'pmdi',
                  'anchor': 'ether_C2',
                  'cap': 'ethyl',
                  'ring': 'H',
                  'template_atom': None},
  'comparison_contrast': 'linker_C2',
  'hypothesis': 'Measure linker_C2 on pmdi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_ndi_dd8125d2de65',
  'name': 'ndi / direct_N / npropyl / H',
  'family': 'imide (n-type)',
  'frag': 'CCCN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N([*:1])C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_ndi_ethyl',
  'search_axes': {'core': 'ndi',
                  'anchor': 'direct_N',
                  'cap': 'npropyl',
                  'ring': 'H',
                  'template_atom': None},
  'comparison_contrast': 'cap_length',
  'hypothesis': 'Measure cap_length on ndi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_ndi_d2ac23774405',
  'name': 'ndi / direct_N / ethyl / F',
  'family': 'imide (n-type)',
  'frag': 'CCN1C(=O)c2ccc3c4c(c(F)cc(c24)C1=O)C(=O)N([*:1])C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_ndi_ethyl',
  'search_axes': {'core': 'ndi',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'F',
                  'template_atom': 5},
  'comparison_contrast': 'ring_F',
  'hypothesis': 'Measure ring_F on ndi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_ndi_da4d3784d1e2',
  'name': 'ndi / direct_N / ethyl / CN',
  'family': 'imide (n-type)',
  'frag': 'CCN1C(=O)c2ccc3c4c(c(C#N)cc(c24)C1=O)C(=O)N([*:1])C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_ndi_ethyl',
  'search_axes': {'core': 'ndi',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'CN',
                  'template_atom': 5},
  'comparison_contrast': 'ring_CN',
  'hypothesis': 'Measure ring_CN on ndi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_ndi_da0ff84e889c',
  'name': 'ndi / ether_C2 / ethyl / H',
  'family': 'imide (n-type)',
  'frag': 'CCN1C(=O)c2ccc3c4c(ccc(c24)C1=O)C(=O)N(CCO[*:1])C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_ndi_ethyl',
  'search_axes': {'core': 'ndi',
                  'anchor': 'ether_C2',
                  'cap': 'ethyl',
                  'ring': 'H',
                  'template_atom': None},
  'comparison_contrast': 'linker_C2',
  'hypothesis': 'Measure linker_C2 on ndi against the unchanged ethyl graft; no performance '
                'improvement assumed.'}]


# Matched-comparison hypotheses registered 2026-10-18; exact-graft properties remain pending.
GROUPS += [{'id': 'search_viologen_1ad45438b9e4',
  'name': 'viologen / direct_N / ethyl / OMe',
  'family': 'pyridine-multi-e',
  'frag': 'CC[n+]1ccc(-c2cc[n+]([*:1])c(OC)c2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'ethylviologen',
  'search_axes': {'core': 'viologen',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'OMe',
                  'template_atom': 2},
  'comparison_contrast': 'ring_OMe',
  'hypothesis': 'Measure ring_OMe on viologen against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_viologen_e1199c9f4288',
  'name': 'viologen / direct_N / ethyl / Me',
  'family': 'pyridine-multi-e',
  'frag': 'CC[n+]1ccc(-c2cc[n+]([*:1])c(C)c2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'ethylviologen',
  'search_axes': {'core': 'viologen',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'Me',
                  'template_atom': 2},
  'comparison_contrast': 'ring_Me',
  'hypothesis': 'Measure ring_Me on viologen against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_pmdi_6e49edb06a78',
  'name': 'pmdi / direct_N / ethyl / OMe',
  'family': 'imide (n-type)',
  'frag': 'CCn1c(=O)c2cc3c(=O)n([*:1])c(=O)c3c(OC)c2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_pmdi_ethyl',
  'search_axes': {'core': 'pmdi',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'OMe',
                  'template_atom': 5},
  'comparison_contrast': 'ring_OMe',
  'hypothesis': 'Measure ring_OMe on pmdi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_pmdi_6c8c35f7efd1',
  'name': 'pmdi / direct_N / ethyl / Me',
  'family': 'imide (n-type)',
  'frag': 'CCn1c(=O)c2cc3c(=O)n([*:1])c(=O)c3c(C)c2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_pmdi_ethyl',
  'search_axes': {'core': 'pmdi',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'Me',
                  'template_atom': 5},
  'comparison_contrast': 'ring_Me',
  'hypothesis': 'Measure ring_Me on pmdi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_ndi_9b5f09730b83',
  'name': 'ndi / direct_N / ethyl / OMe',
  'family': 'imide (n-type)',
  'frag': 'CCN1C(=O)c2ccc3c4c(c(OC)cc(c24)C1=O)C(=O)N([*:1])C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_ndi_ethyl',
  'search_axes': {'core': 'ndi',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'OMe',
                  'template_atom': 5},
  'comparison_contrast': 'ring_OMe',
  'hypothesis': 'Measure ring_OMe on ndi against the unchanged ethyl graft; no performance '
                'improvement assumed.'},
 {'id': 'search_ndi_573f69a7895b',
  'name': 'ndi / direct_N / ethyl / Me',
  'family': 'imide (n-type)',
  'frag': 'CCN1C(=O)c2ccc3c4c(c(C)cc(c24)C1=O)C(=O)N([*:1])C3=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'family-balanced computational comparison',
                 'exact_grafted_experiment': False,
                 'production_qm': False,
                 'route_status': 'Exact precursor and graft route unverified. Secondary-alkyl '
                                 'viologen caps and mono-N-H imide preparation require particular '
                                 'route review.',
                 'source_pool': '/nfs/lambda_stor_01/homes/rzhu/0_redox/results/discovery/systematic_families_2026-10-18/structures.json',
                 'source_pool_sha256': '17c8a9acbe62a81b902c8cc1bed443c43fec22f85d4e391b26b6f3a6952545ba',
                 'parent_metrics_inherited': False},
  'comparison_baseline': 'design_ndi_ethyl',
  'search_axes': {'core': 'ndi',
                  'anchor': 'direct_N',
                  'cap': 'ethyl',
                  'ring': 'Me',
                  'template_atom': 5},
  'comparison_contrast': 'ring_Me',
  'hypothesis': 'Measure ring_Me on ndi against the unchanged ethyl graft; no performance '
                'improvement assumed.'}]


# Matched-comparison hypotheses registered 2026-10-19; exact-graft properties remain pending.
GROUPS += [{'id': 'hyp_cata_ndi_ethyl',
  'name': 'N-(4-methylbenzyl)-N-ethyl 1,2,3,4-naphthalene diimide',
  'family': 'imide (n-type)',
  'frag': 'CCn1c(=O)c2c3ccccc3c3c(=O)n([*:1])c(=O)c3c2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'hypothesis': 'Test cata vs conventional NDI as a core-topology isomer pair with the same '
                'formula, ethyl cap and N-graft anchor.',
  'comparison_baseline': 'design_ndi_ethyl',
  'provenance': {'type': 'source-directed alternative imide core',
                 'publication': 'https://doi.org/10.3762/bjoc.20.155',
                 'source_parent': '7-Ph',
                 'exact_grafted_experiment': False,
                 'parent_metrics_inherited': False,
                 'production_qm': False,
                 'route_status': 'Source parents are bis-phenyl/bis-hexyl. Mixed N-ethyl/N-graft '
                                 'precursor isolation and selectivity are unverified.',
                 'registered': '2026-10-19',
                 'source': 'https://doi.org/10.3762/bjoc.20.155'},
  'exact_grafted_smiles': 'CCn1c(=O)c2c3ccccc3c3c(=O)n(Cc4ccc(C)cc4)c(=O)c3c2c1=O',
  'inchikey': 'NANYKQXYZMYUQE-UHFFFAOYSA-N',
  'formula': 'C24H18N2O4',
  'total_atoms_including_H': 48,
  'MW': 398.41800000000023,
  'precursor_smiles': 'CCn1c(=O)c2c3ccccc3c3c(=O)[nH]c(=O)c3c2c1=O',
  'chemistry_notes': 'Test cata vs conventional NDI as a core-topology isomer pair with the same '
                     'formula, ethyl cap and N-graft anchor. Source parents are '
                     'bis-phenyl/bis-hexyl. Mixed N-ethyl/N-graft precursor isolation and '
                     'selectivity are unverified. Parent measurements are in CH2Cl2 and are not '
                     'exact-graft or MeCN labels.'},
 {'id': 'hyp_cata_adi_ethyl',
  'name': 'N-(4-methylbenzyl)-N-ethyl 1,2,3,4-anthracene diimide',
  'family': 'imide (n-type)',
  'frag': 'CCn1c(=O)c2c3cc4ccccc4cc3c3c(=O)n([*:1])c(=O)c3c2c1=O',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'hypothesis': 'Test adding one fused benzene to the cata-NDI core, with the same ethyl cap and '
                'N-graft anchor; core size changes explicitly.',
  'comparison_baseline': 'hyp_cata_ndi_ethyl',
  'provenance': {'type': 'source-directed alternative imide core',
                 'publication': 'https://doi.org/10.3762/bjoc.20.155',
                 'source_parent': '8-Ph',
                 'exact_grafted_experiment': False,
                 'parent_metrics_inherited': False,
                 'production_qm': False,
                 'route_status': 'Source parents are bis-phenyl/bis-hexyl. Mixed N-ethyl/N-graft '
                                 'precursor isolation and selectivity are unverified.',
                 'registered': '2026-10-19',
                 'source': 'https://doi.org/10.3762/bjoc.20.155'},
  'exact_grafted_smiles': 'CCn1c(=O)c2c3cc4ccccc4cc3c3c(=O)n(Cc4ccc(C)cc4)c(=O)c3c2c1=O',
  'inchikey': 'LKTOOBTUSSFQCC-UHFFFAOYSA-N',
  'formula': 'C28H20N2O4',
  'total_atoms_including_H': 54,
  'MW': 448.47800000000024,
  'precursor_smiles': 'CCn1c(=O)c2c3cc4ccccc4cc3c3c(=O)[nH]c(=O)c3c2c1=O',
  'chemistry_notes': 'Test adding one fused benzene to the cata-NDI core, with the same ethyl cap '
                     'and N-graft anchor; core size changes explicitly. Source parents are '
                     'bis-phenyl/bis-hexyl. Mixed N-ethyl/N-graft precursor isolation and '
                     'selectivity are unverified. Parent measurements are in CH2Cl2 and are not '
                     'exact-graft or MeCN labels.'}]


# Matched-comparison hypotheses registered 2026-10-19; exact-graft properties remain pending.
GROUPS += [{'id': 'study_aryl_0d7aa75989e9',
  'name': 'AV-3a ester at aryl C3',
  'family': 'pyridine-multi-e',
  'frag': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O[*:1])c4)cc3)cc2)c1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3a',
                 'source_SI_page': 9,
                 'source_parent_component_smiles': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O)c4)cc3)cc2)c1',
                 'source_parent_salt_smiles': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O)c4)cc3)cc2)c1.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [3]}},
 {'id': 'study_aryl_45a6c406cdc7',
  'name': 'AV-3b ester at aryl C4',
  'family': 'pyridine-multi-e',
  'frag': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O[*:1])cc4)cc3)cc2)cc1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3b',
                 'source_SI_page': 10,
                 'source_parent_component_smiles': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O)cc4)cc3)cc2)cc1',
                 'source_parent_salt_smiles': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O)cc4)cc3)cc2)cc1.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [4]}},
 {'id': 'study_aryl_c285539322fa',
  'name': 'AV-3c ester at aryl C3/5',
  'family': 'pyridine-multi-e',
  'frag': 'O=C(O)c1cc(C(=O)O)cc(-[n+]2ccc(-c3cc[n+](-c4cc(C(=O)O)cc(C(=O)O[*:1])c4)cc3)cc2)c1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3c',
                 'source_SI_page': 10,
                 'source_parent_component_smiles': 'O=C(O)c1cc(C(=O)O)cc(-[n+]2ccc(-c3cc[n+](-c4cc(C(=O)O)cc(C(=O)O)c4)cc3)cc2)c1',
                 'source_parent_salt_smiles': 'O=C(O)c1cc(C(=O)O)cc(-[n+]2ccc(-c3cc[n+](-c4cc(C(=O)O)cc(C(=O)O)c4)cc3)cc2)c1.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [3, 5]}},
 {'id': 'study_aryl_4602f4a856c7',
  'name': 'AV-3d ester at aryl C3',
  'family': 'pyridine-multi-e',
  'frag': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O)c(C(=O)O[*:1])c4)cc3)cc2)cc1C(=O)O',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3d',
                 'source_SI_page': 11,
                 'source_parent_component_smiles': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O)c(C(=O)O)c4)cc3)cc2)cc1C(=O)O',
                 'source_parent_salt_smiles': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O)c(C(=O)O)c4)cc3)cc2)cc1C(=O)O.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [3]}},
 {'id': 'study_aryl_e2795f06393d',
  'name': 'AV-3d ester at aryl C4',
  'family': 'pyridine-multi-e',
  'frag': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O[*:1])c(C(=O)O)c4)cc3)cc2)cc1C(=O)O',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3d',
                 'source_SI_page': 11,
                 'source_parent_component_smiles': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O)c(C(=O)O)c4)cc3)cc2)cc1C(=O)O',
                 'source_parent_salt_smiles': 'O=C(O)c1ccc(-[n+]2ccc(-c3cc[n+](-c4ccc(C(=O)O)c(C(=O)O)c4)cc3)cc2)cc1C(=O)O.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [4]}},
 {'id': 'study_aryl_c88a8d5645d4',
  'name': 'AV-3e ester at aryl C2',
  'family': 'pyridine-multi-e',
  'frag': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O)c4C(=O)O[*:1])cc3)cc2)c1C(=O)O',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3e',
                 'source_SI_page': 12,
                 'source_parent_component_smiles': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O)c4C(=O)O)cc3)cc2)c1C(=O)O',
                 'source_parent_salt_smiles': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O)c4C(=O)O)cc3)cc2)c1C(=O)O.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [2]}},
 {'id': 'study_aryl_693080e8c287',
  'name': 'AV-3e ester at aryl C3',
  'family': 'pyridine-multi-e',
  'frag': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O[*:1])c4C(=O)O)cc3)cc2)c1C(=O)O',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3e',
                 'source_SI_page': 12,
                 'source_parent_component_smiles': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O)c4C(=O)O)cc3)cc2)c1C(=O)O',
                 'source_parent_salt_smiles': 'O=C(O)c1cccc(-[n+]2ccc(-c3cc[n+](-c4cccc(C(=O)O)c4C(=O)O)cc3)cc2)c1C(=O)O.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [3]}},
 {'id': 'study_aryl_9992bd3479e1',
  'name': 'AV-3f ether at aryl C2',
  'family': 'pyridine-multi-e',
  'frag': 'Oc1ccc(O)c(-[n+]2ccc(-c3cc[n+](-c4cc(O)ccc4O[*:1])cc3)cc2)c1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3f',
                 'source_SI_page': 15,
                 'source_parent_component_smiles': 'Oc1ccc(O)c(-[n+]2ccc(-c3cc[n+](-c4cc(O)ccc4O)cc3)cc2)c1',
                 'source_parent_salt_smiles': 'Oc1ccc(O)c(-[n+]2ccc(-c3cc[n+](-c4cc(O)ccc4O)cc3)cc2)c1.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [2]}},
 {'id': 'study_aryl_6f1ff24b8d40',
  'name': 'AV-3f ether at aryl C5',
  'family': 'pyridine-multi-e',
  'frag': 'Oc1ccc(O)c(-[n+]2ccc(-c3cc[n+](-c4cc(O[*:1])ccc4O)cc3)cc2)c1',
  'states': [['ox2', 2, 1, 0], ['ox1', 1, 2, -1], ['neu', 0, 1, -2]],
  'hypothesis': 'Retain both N-aryl bonds while varying O-attachment; exact-graft behavior '
                'unverified.',
  'provenance': {'type': 'literature-parent-derived computational hypothesis',
                 'source_doi': '10.1021/acsaem.5c02255',
                 'source_parent_id': 'AV-3f',
                 'source_SI_page': 15,
                 'source_parent_component_smiles': 'Oc1ccc(O)c(-[n+]2ccc(-c3cc[n+](-c4cc(O)ccc4O)cc3)cc2)c1',
                 'source_parent_salt_smiles': 'Oc1ccc(O)c(-[n+]2ccc(-c3cc[n+](-c4cc(O)ccc4O)cc3)cc2)c1.[Cl-].[Cl-]',
                 'parent_metrics_inherited': False,
                 'exact_grafted_experiment': False,
                 'protonation_model': 'As-isolated acid/phenol component, modeled as bare dication '
                                      'in MeCN; aqueous pH14 speciation not inferred',
                 'route_status': 'Parent synthesis reported; exact Merrifield graft synthesis and '
                                 'its stability unverified',
                 'source_attachment_positions': [5]}}]


# New two-electron scaffold families (expansion 2026-11-16). The declared
# viologen/PMDI/NDI substituent space was saturated, so these add validated 2e CORES
# beyond it: 5,10-dihydrophenazine (a 2e p-type CATHOLYTE, the side the pool lacked)
# and 9,10-anthraquinone (a 2e n-type anolyte). Grafted via the standard 4-methylbenzyl
# tether proxy. Exact graft properties are uncomputed (production_qm False); these are
# definitions for the route-review + QM queue, not inherited parent measurements.
GROUPS += [
 {'id': 'exp_dmpz_mb_ethyl',
  'name': "N-(4-methylbenzyl)-N'-ethyl-5,10-dihydrophenazine",
  'family': 'phenazine (p-type)',
  'frag': '[*:1]N1c2ccccc2N(CC)c2ccccc21',
  'states': [['neu', 0, 1, 0], ['ox1', 1, 2, 1], ['ox2', 2, 1, 2]],
  'provenance': {'type': 'new 2e scaffold family (expansion 2026-11-16)',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'exact_grafted_experiment': False,
                 'direct_merrifield_route_verified': False,
                 'parent_metrics_inherited': False,
                 'source': 'https://www.sciencedirect.com/science/article/pii/S2451929419303158',
                 'family_grounding': 'Parent N,N-dialkyl-5,10-dihydrophenazine (BMEPZ) is a '
                                     'two-electron nonaqueous catholyte: two 1e oxidations at '
                                     '-0.29 and +0.50 V vs Fc in MeCN, 0.5 M solubility, 99.94%/'
                                     'cycle over 200 cycles (Bio-inspired Molecular Redesign, 2019).',
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Fills the missing 2e p-type catholyte slot. N-H alkylation on 5,10-'
                     'dihydrophenazine; ethyl is the matched control cap. Second oxidation is '
                     'validated for the parent family (unlike the removed phenothiazine 2e). '
                     'Exact graft synthesis and its stability require separate route review.',
  'review_priority': 'catholyte 2e lead (new family)'},
 {'id': 'hyp_dmpz_mb_methoxyethyl',
  'name': "N-(4-methylbenzyl)-N'-(2-methoxyethyl)-5,10-dihydrophenazine",
  'family': 'phenazine (p-type)',
  'frag': '[*:1]N1c2ccccc2N(CCOC)c2ccccc21',
  'states': [['neu', 0, 1, 0], ['ox1', 1, 2, 1], ['ox2', 2, 1, 2]],
  'comparison_baseline': 'exp_dmpz_mb_ethyl',
  'hypothesis': 'Methoxyethyl cap mirrors the validated BMEPZ solubilizing substituent; test '
                'solubility/mass against the ethyl control. No performance improvement assumed.',
  'provenance': {'type': 'new 2e scaffold family (expansion 2026-11-16)',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'exact_grafted_experiment': False,
                 'direct_merrifield_route_verified': False,
                 'parent_metrics_inherited': False,
                 'source': 'https://www.sciencedirect.com/science/article/pii/S2451929419303158',
                 'family_grounding': 'BMEPZ (the 2-methoxyethyl-capped parent) is the measured 2e '
                                     'catholyte; this graft substitutes one cap with the Merrifield tether.',
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Dihydrophenazine 2e catholyte with the BMEPZ-type glycol cap for MeCN '
                     'solubility. Exact graft synthesis and stability require separate route review.',
  'review_priority': 'catholyte 2e (new family)'},
 {'id': 'hyp_aq_2_benzyloxy',
  'name': '2-[(4-methylbenzyl)oxy]-9,10-anthraquinone',
  'family': 'quinone (n-type)',
  'frag': '[*:1]Oc1ccc2C(=O)c3ccccc3C(=O)c2c1',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'provenance': {'type': 'new 2e scaffold family (expansion 2026-11-16)',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'exact_grafted_experiment': False,
                 'direct_merrifield_route_verified': False,
                 'parent_metrics_inherited': False,
                 'family_grounding': '9,10-anthraquinone is a tier-A two-wave 2e anolyte in the '
                                     'in-repo MeCN benchmark (config/benchmark_mecn.py).',
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Anthraquinone grafted through a 2-position ether (benzylic-Cl displacement by '
                     '2-hydroxyanthraquinone). Extends the quinone 2e anolyte family to the larger '
                     'anthraquinone core. Exact graft synthesis and stability require route review.',
  'review_priority': 'anolyte 2e (new core)'},
 {'id': 'hyp_aq_2_amino',
  'name': '2-[(4-methylbenzyl)amino]-9,10-anthraquinone',
  'family': 'quinone (n-type)',
  'frag': '[*:1]Nc1ccc2C(=O)c3ccccc3C(=O)c2c1',
  'states': [['neu', 0, 1, 0], ['red1', -1, 2, -1], ['red2', -2, 1, -2]],
  'comparison_baseline': 'hyp_aq_2_benzyloxy',
  'hypothesis': 'Amino vs ether anchor on the anthraquinone 2-position; test the donor effect on '
                'E1/E2 and the two-wave spacing. No performance improvement assumed.',
  'provenance': {'type': 'new 2e scaffold family (expansion 2026-11-16)',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'exact_grafted_experiment': False,
                 'direct_merrifield_route_verified': False,
                 'parent_metrics_inherited': False,
                 'family_grounding': '9,10-anthraquinone is a tier-A two-wave 2e anolyte in the '
                                     'in-repo MeCN benchmark (config/benchmark_mecn.py).',
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Anthraquinone grafted through a 2-position secondary amine (benzylic-Cl '
                     'displacement by 2-aminoanthraquinone). Exact graft synthesis and stability '
                     'require route review.',
  'review_priority': 'anolyte 2e (new core)'}]

# Catholyte-focused expansion (2026-11-17). The full-pool scorecard shows the catholyte
# Pareto front has only 3 members vs 7 anolytes, and all 3 are dihydrophenazine/phenothiazine
# — the 2e catholyte side is the real scarcity. These deepen the ONE validated 2e catholyte
# family (5,10-dihydrophenazine, BMEPZ parent) along two grounded levers: (a) a smaller N'-cap
# raises gravimetric capacity (less dead mass), (b) N'/ring electron-withdrawing vs -donating
# substituents tune the oxidation potential (EWG raise E -> higher cell voltage; EDG lower E,
# usually better radical-cation stability). Same 4-methylbenzyl Merrifield tether; properties
# uncomputed (production_qm False) — definitions for the QM queue, not inherited measurements.
GROUPS += [
 {'id': 'exp_dmpz_mb_methyl',
  'name': "N-(4-methylbenzyl)-N'-methyl-5,10-dihydrophenazine",
  'family': 'phenazine (p-type)',
  'frag': '[*:1]N1c2ccccc2N(C)c2ccccc21',
  'states': [['neu', 0, 1, 0], ['ox1', 1, 2, 1], ['ox2', 2, 1, 2]],
  'comparison_baseline': 'exp_dmpz_mb_ethyl',
  'hypothesis': 'Smallest alkyl cap on the validated DHP core -> lowest dead mass -> highest '
                'gravimetric capacity of the DHP catholyte series (~91 vs 89 mAh/g for ethyl). '
                'No change to the redox chemistry assumed.',
  'provenance': {'type': 'catholyte 2e expansion (2026-11-17)',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'exact_grafted_experiment': False,
                 'direct_merrifield_route_verified': False,
                 'parent_metrics_inherited': False,
                 'source': 'https://www.sciencedirect.com/science/article/pii/S2451929419303158',
                 'family_grounding': 'N,N-dialkyl-5,10-dihydrophenazine (BMEPZ) is a measured 2e '
                                     'nonaqueous catholyte (-0.29/+0.50 V vs Fc, MeCN).',
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Capacity-maximizing member of the DHP 2e catholyte family. N-H methylation '
                     'on 5,10-dihydrophenazine. Route/stability require separate review.',
  'review_priority': 'catholyte 2e (capacity lead)'},
 {'id': 'hyp_dmpz_mb_trifluoroethyl',
  'name': "N-(4-methylbenzyl)-N'-(2,2,2-trifluoroethyl)-5,10-dihydrophenazine",
  'family': 'phenazine (p-type)',
  'frag': '[*:1]N1c2ccccc2N(CC(F)(F)F)c2ccccc21',
  'states': [['neu', 0, 1, 0], ['ox1', 1, 2, 1], ['ox2', 2, 1, 2]],
  'comparison_baseline': 'exp_dmpz_mb_ethyl',
  'hypothesis': 'Electron-withdrawing 2,2,2-trifluoroethyl cap lowers HOMO density -> raises both '
                'oxidation potentials vs the ethyl control -> higher full-cell voltage. Tests '
                'whether the 2nd oxidation stays reversible under an EWG. No value assumed.',
  'provenance': {'type': 'catholyte 2e expansion (2026-11-17)',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'exact_grafted_experiment': False,
                 'direct_merrifield_route_verified': False,
                 'parent_metrics_inherited': False,
                 'source': 'https://www.sciencedirect.com/science/article/pii/S2451929419303158',
                 'family_grounding': 'BMEPZ 2e catholyte parent; substituent-potential trend '
                                     '(EWG raises oxidation potential) is standard for '
                                     'N-aryl/N-alkyl dihydrophenazine redox.',
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Potential-raising member of the DHP catholyte family via an N-CH2CF3 cap. '
                     'Route/stability require separate review.',
  'review_priority': 'catholyte 2e (voltage lead)'},
 {'id': 'hyp_dmpz_27dimethyl_ethyl',
  'name': "N-(4-methylbenzyl)-N'-ethyl-2,7-dimethyl-5,10-dihydrophenazine",
  'family': 'phenazine (p-type)',
  'frag': '[*:1]N1c2cc(C)ccc2N(CC)c2ccc(C)cc21',
  'states': [['neu', 0, 1, 0], ['ox1', 1, 2, 1], ['ox2', 2, 1, 2]],
  'comparison_baseline': 'exp_dmpz_mb_ethyl',
  'hypothesis': 'Ring methyl donors at 2,7 lower the oxidation potentials modestly and can '
                'stabilize the oxidized states sterically/electronically vs the unsubstituted '
                'ethyl control. Tests the ring-donor lever. No value assumed.',
  'provenance': {'type': 'catholyte 2e expansion (2026-11-17)',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'exact_grafted_experiment': False,
                 'direct_merrifield_route_verified': False,
                 'parent_metrics_inherited': False,
                 'source': 'https://www.sciencedirect.com/science/article/pii/S2451929419303158',
                 'family_grounding': 'BMEPZ 2e catholyte parent; ring-EDG lowers oxidation '
                                     'potential in dihydrophenazine/phenazine redox.',
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Ring-substituted DHP (2,7-dimethyl) probing the core-donor effect on E1/E2. '
                     'Route/stability require separate review.',
  'review_priority': 'catholyte 2e (ring-donor probe)'},
 {'id': 'hyp_dmpz_27dimethoxy_ethyl',
  'name': "N-(4-methylbenzyl)-N'-ethyl-2,7-dimethoxy-5,10-dihydrophenazine",
  'family': 'phenazine (p-type)',
  'frag': '[*:1]N1c2cc(OC)ccc2N(CC)c2ccc(OC)cc21',
  'states': [['neu', 0, 1, 0], ['ox1', 1, 2, 1], ['ox2', 2, 1, 2]],
  'comparison_baseline': 'hyp_dmpz_27dimethyl_ethyl',
  'hypothesis': 'Stronger 2,7-methoxy donors lower the oxidation potentials further and typically '
                'improve radical-cation stability (as in dialkoxy p-type catholytes). Tests the '
                'strong ring-donor limit. No value assumed.',
  'provenance': {'type': 'catholyte 2e expansion (2026-11-17)',
                 'measured_for_this_structure': False,
                 'production_qm': False,
                 'exact_grafted_experiment': False,
                 'direct_merrifield_route_verified': False,
                 'parent_metrics_inherited': False,
                 'source': 'https://www.sciencedirect.com/science/article/pii/S2451929419303158',
                 'family_grounding': 'BMEPZ 2e catholyte parent; methoxy ring-donors lower '
                                     'oxidation potential and stabilize radical cations in '
                                     'dihydrophenazine/dialkoxyarene p-type redox.',
                 'registration_stage': 'exact structure/state definition; properties pending'},
  'chemistry_notes': 'Ring-substituted DHP (2,7-dimethoxy) probing the strong core-donor limit. '
                     'Route/stability require separate review.',
  'review_priority': 'catholyte 2e (ring-donor probe)'}]
