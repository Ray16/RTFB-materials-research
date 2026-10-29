# Preliminary DBUH–chloride association check

[Reproducible comparison](comparison.json) · [Calculation](reproduce.py)

Goodwin et al. report a DBUH+–chloride association constant of **419 M⁻¹** in
[their article](https://doi.org/10.1039/D1OB00282A), note 48. The downloaded
[author-accepted manuscript](../../calcs/validation/dbu_chloride_anchor_2026-10-29/goodwin_accepted_manuscript.pdf)
(SHA-256 `71510f845d5f6f050f9aeff8e5c8f918b19ca3acf5a77f12463060e164bff4ff`)
confirms the value and directs readers to its supporting information. The
publisher's SI URL returns HTTP 404 at this check. Search-indexed SI text states
CD3CN and DBUH·PF6; it references figures S66/S67 and a general titration
protocol. **Temperature, binding-fit details, experimental uncertainty, water
content and counterion/activity corrections for this exact measurement remain
unverified.** The 298 K caption for a different receptor is not evidence of the
DBUH experiment's temperature.

A conditional conversion at **assumed 298.15 K**, with a 1 M standard
concentration, yields −RT ln(419 M⁻¹ × 1 M) = **−14.97 kJ/mol**. The archived
base-assisted calculation reports **−4.80 kJ/mol** for one DBUH+/chloride
contact pair relative to separated ions, corresponding formally to about
6.9 M⁻¹ at the same assumed temperature and standard state. The difference is
approximately **+10.17 kJ/mol, modeled minus conditional experimental value**.
The producer's stored contact and separated corrections are checked to differ by
the reported association term. This does not independently rerun or physically
revalidate its underlying quantum states.

This is a preliminary diagnostic, not an error bar, calibration, or proof of
incorrect grafting free energies. Experiment can include solution populations
and competition with PF6− and other ions; the calculation represents one
optimized contact-pair state, and CD3CN is not exactly modeled CH3CN. Ionic
activities, finite concentration, solvent water, temperature and conformer
sampling could affect comparison. The association shift is shared under a fixed
DBU convention, so it cannot change within-base candidate ordering, and it
cannot be blindly added to the full reaction: there may be error cancellation
and additional solution species.

The existing imide and phenolic grafting cases have **no identified matched
full-reaction equilibrium measurement** in this focused search. An [imide
preparation](https://www.orgsyn.org/demo.aspx?prep=CV2P0083) establishes a
synthetic route under different conditions; isolated yield is not equilibrium
ΔG. [Phenol/phenolate solution equilibria](https://doi.org/10.1016/S0021-9614(02)00323-3)
illustrate competing homoconjugation, but are not measurements of our grafting
reactions. The [earlier acid-capture comparison](../../../archive/0_redox/results/validation/acid_base_external_2026-10-22/README.md)
and [DIPEA comparison](../../../archive/0_redox/results/validation/dipea_uv_anchor_2026-10-26/README.md)
remain separate controls. No candidate admission, full-reaction benchmark or
new DFT calculation is claimed.

Next: obtain the SI and inspect the exact fit, temperature, concentrations and
salt conventions; seek independently measured imide/phenol full-reaction or
matched component equilibria. Once the source convention is resolved, compare
multiple DBUH chloride arrangements and relevant solvent/counterion models with
unchanged physical QC. Do not infer experimental grafting yield from ΔG.
