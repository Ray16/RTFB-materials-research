# External check of a defined acid-capture term

The [comparison](comparison.json) checks `B + HCl → BH+ + Cl−` in MeCN against
author-tabulated pKa anchors. It is a scoped, preliminary thermodynamic check;
source temperature/activity details still require review. It does not establish
the accuracy of a complete grafting reaction.

| DBU separated-ion capture | kJ/mol |
| --- | ---: |
| All-DFT term | −102.78 |
| Term from pKaH(DBU) = 24.31 and pKa(HCl) = 10.30 at 298.15 K | −79.97 |
| DFT minus source-derived term | −22.81 |

The base is matched by its exact neutral molecular graph to the
[author-maintained basicity collection](https://analytical.chem.ut.ee/databases/pka-values/),
which uses the revised MeCN scale described in
[Tshepelevitsh et al.](https://doi.org/10.1002/ejoc.201900956).
The [author acid table](https://analytical.chem.ut.ee/HA_UT/pKa_Table_Water_Acetonitrile_DCE_Gas_Phase.pdf)
places HCl at 10.30 in the MeCN column and links
[Raamat et al.](https://doi.org/10.1002/poc.2946) for qualifications.
That PDF row and column were checked visually. This is an author-recommended
anchor, not an assertion of a directly matched experimental condition.
Raw downloads and hashes are retained under
`calcs/validation/acid_base_external_2026-10-22/sources/`.

The pKa conversion uses `ΔG = −RT ln(10) [pKa(BH+) − pKa(HCl)]`.
The number of separately translating solutes is unchanged, so the 1-atm/1-M
reaction correction cancels. No isolated-proton free energy is introduced.

Replacing only this term while keeping the DFT HCl-reference reaction gives BQ
values of −22.33 and −23.94 kJ/mol with DBU, compared with the original all-DFT
−45.14 and −46.75 kJ/mol. Both signs remain favorable in that limited sensitivity
calculation. These are not calibrated predictions: HCl cancels in the complete
all-DFT cycle, so term errors and shared-error cancellation cannot be assigned to
individual species or transferred uncritically to the whole reaction.

No exact DIPEA graph was found in the checked collection. Its external comparison
remains blank. Diisopropylamine and triethylamine are distinct molecules and are
not substitutes. No contact-pair association value is inferred from free-ion pKa.
Existing computed grafting values and synthesis-route holds remain unchanged.
