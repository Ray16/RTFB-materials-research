# Vinyl monomers for functionalization before polymerization

The confirmed synthesis route is **attach the redox center to a vinyl monomer,
then polymerize**. These isolated proposals replace the local model's para methyl
cap with a terminal vinyl group. The target scaffold is pure para-substituted
4-vinylbenzyl; other positional isomers require separate proposals.

The [proposal records](proposals.json) and [CSV](proposals.csv) cover all 75 design
graphs and the ethyl-viologen comparison control: 76 unique vinyl-bearing graphs.
Every graph passes the existing attachment, stereochemistry, charge-ladder and
electron-parity contracts. A separate [review](independent_review.json) reconstructs
the original capped graphs through an independent reaction and checks elemental
mass and nominal-capacity arithmetic.

These are graph proposals, **not computed vinyl-monomer properties or verified
synthesis routes**. Redox potentials and solution reorganization energies are
blank. Capped-model values are not transferred. The declared charge ladder is a
proposal to test, not a spin or thermodynamic result for the vinyl target.
No shared manifest, configuration, conformer or canonical calculation is changed.

The [isolated initial pilot](pilot/preparation.json) contains three explicit
vinyl targets and nine declared state/seed sets. Its generated conformers pass
81 geometry/hydrogen-topology seed checks. These are starting inputs; neither
ML relaxation nor solution quantum properties are implied by that preparation.

The [bounded solution pilot](pilot/README.md) has now started. Six actual DFT
workers were independently verified on reserved H100 GPUs; nine charge-state
calculations are queued. Results remain unpublished until state, spin and
complete-couple checks admit them. See the timestamped launch and observation
records for the recorded startup state.

The initial modeling and route-review controls are the vinyl counterparts of
conventional NDI ethyl, PMDI ethyl and ethyl viologen. They provide established
core families against which more exploratory cata-imides and substituted viologens
can be compared. This pilot selection is not a claim of experimental superiority.

The next checks are exact precursor/attachment chemistry, preservation of the
vinyl group, polymerization compatibility, and redox-center retention. Quinone
branches require particular review of radical trapping/inhibition and O/C
attachment selectivity. Relevant precedents and limitations are linked in the
[maintained figure guide](../../candidates/figures/README.md).

Molecular masses and nominal capacities refer to the charge-bearing component
where ionic. Counterion choice, comonomer fraction, polymer loading, swelling,
charge transport and electrochemical accessibility remain unqualified. The
cationic viologen graph does not represent an isolated counterion-free substance.
Polymerization changes the electronic environment; neither the capped model nor
the isolated vinyl monomer alone establishes polymer redox performance.

Preparation: `python scripts/pipeline/prepare_vinyl_monomers.py --out <new-output-folder>`.
Review this prepared set: `python results/discovery/vinyl_monomers/review_preparation.py`.
The preparer refuses to overwrite an existing proposal folder. It submits no jobs.
