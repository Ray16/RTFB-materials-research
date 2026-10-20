# Practical attachment of redox motifs to Merrifield polymers

Reviewed 19 October 2026. Scope: functionalize an existing chloromethylated
polystyrene, including crosslinked Merrifield beads. This is a focused literature
assessment and modeling guide, not a validated experimental protocol for our
75 exact structures.

## What the attachment reaction actually is

The polymer supplies a benzyl chloride handle, `P–CH2–Cl`. A nucleophile replaces
chloride at that benzyl carbon. The redox motif therefore needs an accessible,
selectively reactive attachment group; a favorable redox potential does not
provide one. Preparing the nucleophile, getting it into the swollen polymer,
and preserving the redox center all belong in the route assessment.

The distinction between neutral and activated precursors matters:

```text
Preformed salt:    P–CH2–Cl + K[Nu]  → P–CH2–Nu + KCl
Base-assisted:    P–CH2–Cl + H–Nu + B → P–CH2–Nu + [BH]+ + Cl−
Pyridyl N:       P–CH2–Cl + :N–R     → [P–CH2–N–R]+ + Cl−
```

These are bookkeeping schemes. Charges on the rest of a motif and its original
counterions must also be retained. Salt solubility, ion pairing and phases need
an explicit model. A carbonate route needs its own balanced species rather than
an arbitrary replacement of HCl by a constant correction.

## Demonstrated routes and the limits of transfer

| Attachment | Primary evidence and reported conditions | What it establishes for this project |
|---|---|---|
| Viologen, pyridyl N alkylation | Kim, Sanford, Vaid and McNeil (2022): Merrifield beads with 1-methyl-4,4′-bipyridinium iodide, refluxing MeCN, 48 h; subsequent NH4PF6 exchange. [Article](https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.202200149) | A close precedent on crosslinked beads in our selected solvent. Each substituted precursor still needs access/selectivity review; the final counterions differ from the initial chloride/iodide product. |
| Imide N alkylation | Kim, Ju, Park and Kang (2019): potassium phthalimide, 1.5 equivalents per chloromethyl site, DMF, 70 °C, 6 h; approximately complete substitution by proton NMR. [Article](https://pubs.rsc.org/en/content/articlehtml/2019/ra/c9ra01131e) | Strong local chemistry precedent for an imide anion. The support was soluble poly(4-chloromethylstyrene), so this does not establish bead penetration, MeCN performance or the exact NDI/PMDI routes. |
| AQ through a carboxylate-derived ester | Oyaizu, Choi and Nishide (2011): anthraquinone-2-carboxylic acid attached to poly(4-chloromethylstyrene). [Article](https://onlinelibrary.wiley.com/doi/abs/10.1002/pat.1968). Coauthor Choi's [thesis, p. 51](https://waseda.repo.nii.ac.jp/record/19599/files/Honbun-5904.pdf) reports NMP/DBU, 70 °C, 24 h and 80% isolated polymer yield. | A credible alternative AQ attachment design. Isolated mass yield is not site conversion. The benzyl ester has a different molecular graph from our phenolic ether candidates; its properties and linkage stability need separate evaluation. |
| Phenolic O attachment | Huang, Zheng and Long (2010) prepared a Merrifield-supported vanillin derivative before converting its aldehyde to an alcohol. [Publisher abstract](https://link.springer.com/article/10.1007/s12039-010-0023-x) | Supports phenolic attachment as a reaction class. Detailed initial grafting conditions were not verified from accessible publisher full text here. Vanillin is not an exact hydroxy-BQ/AQ precursor. |
| Competing support reaction | Sampson, Simmons and Bradley (2001): HMDS bases could convert Merrifield chloride sites into aminomethyl resin during attempted alkoxide chemistry. [Article](https://www.sciencedirect.com/science/article/abs/pii/S0040403901009947) | Base identity can change the polymer itself. Increasing base strength cannot be treated as a universal grafting remedy. |

The AQ paper reports reversible charge storage in its attached quinone system;
this supports investigating that linker, without establishing our potentials or
long-term operation in the proposed bead-based flow battery.

## What “easy to graft” must include

**Attachment and selectivity.** Review the exact precursor and attachment atom.
An imide route requires the intended free N–H or its salt; fully N-substituted
diimides lack that handle. Multiple reactive groups can produce unintended
attachment patterns or bridges. Pyridyl quaternization uses an available nitrogen
lone pair and does not require an N–H group. Alternative side-chain handles are
possible designs, but changing a linker creates a new candidate graph.

**Polymer accessibility.** Record crosslink content, initial chloride loading,
bead size, reaction solvent and measured swelling. A solution-polymer reaction
does not establish an equally efficient crosslinked-bead reaction. Reaction
solvent and battery electrolyte can be different: the 2022 bead study used DMF
for ferrocene attachment and MeCN for viologen attachment and operation. Its
viologen conversion was estimated from elemental analysis; 75% was an
approximate consistent assumption, not an exact measured grafting yield.

**Useful loading.** Confirm covalent attachment after exhaustive washing,
quantify remaining chloride and retained motif, and test leaching in electrolyte.
Disappearance of one IR band alone is insufficient to assign precise loading.
Use complementary elemental/halogen analysis and appropriate solid-state or
solution spectroscopy. Measure accessible charge by electrochemistry as well.
These are proposed acceptance measurements, not completed experiments here.

For final dry-bead loading `f` in mmol of redox centers per gram, the nominal
capacity is `Q = 26.80 × n × f` mAh/g, where `n` is the electron count used in the
operating window. Use final functionalized mass, including the chosen counterion
convention. Starting chloride loading and molecular capacity do not determine
this value. Swollen bead volume is additionally needed for reservoir energy
density. Our interpretation is that loading, transport and stacking must be
optimized together rather than maximizing graft density alone.

## Design requirements for future molecules

Design each new structure as a precursor-to-product pair, rather than adding an
attachment point after ranking an isolated redox core. Its design record must
include:

1. Exact precursor and grafted-product structures, the intended attachment atom
   and the bond formed at the benzyl-chloride site.
2. A realistic route to the precursor, with competing reactive groups identified
   and a selectivity or protection strategy where necessary.
3. Literature precedent, distinguishing the exact substrate/support from a
   molecular analogue, a family-level analogy or a proposed reaction.
4. Reactive form, solvent, base, counterions, temperature and an atom/charge-
   balanced reaction, with unresolved quantities stated explicitly.
5. Expected compatibility of the redox core and linkage with grafting conditions
   and subsequent charge/discharge; expected polymer accessibility and a plan
   to measure loading.

Our design priority is a selective, accessible handle whose chemistry preserves
the redox core. Linker mass, flexibility, steric shielding and achievable loading
then enter the property comparison. Any proposed anti-stacking substituent must
also be checked for obstruction of the attachment site. Evidence sufficient to
justify computation does not establish an experimentally validated grafting
yield; preserve that distinction in the candidate record.

## Consequences for our computed grafting energies

The current N–H/O–H calculation uses neutral precursor plus capped benzyl
chloride and produces molecular HCl. The positive BQ/AQ values describe that
reference reaction. They do not include base consumption, salt formation or
precipitation, and do not give an activation barrier or grafting yield.

For a specified soluble base, a consistent thermodynamic cycle can add the
free energy of `B + HCl → BH+ + Cl−`. For a preformed potassium nucleophile,
compute the salt reaction with a defined ionic/phase convention instead.
Explicit ion-pair and solvent sensitivity are needed before interpreting small
differences. No new condition-specific free energies were computed in this review.

Keep the highlighted BQ/AQ designs **on hold**, rather than declaring the whole
quinone family ungraftable. The accessible sources reviewed here did not
establish a practical route for those exact phenolic precursors. Prioritize
review of the close viologen precedent, the imide-anion route, and a separately
registered AQ-carboxylate alternative. None is automatically admitted to synthesis
or to a new full-property campaign by this family-level review.

The machine-readable sources are in
[route_evidence.json](../results/grafting/route_evidence.json); candidate-specific
reviews remain separate. The [route screen](../results/grafting/route_screen.json)
requires precursor access, selectivity, defined conditions, a balanced route and
redox-core retention before admission. The
[synthesis-route guide](SYNTHESIS_ROUTE.md) records the prospective evaluation order.
