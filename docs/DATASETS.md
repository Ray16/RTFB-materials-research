# Datasets & Databases of Redox-Active Molecules

Vetted inventory for the MeCN organic-redox screening project (viologen/pyridinium,
TEMPO/nitroxide, phenothiazine, quinone/anthraquinone). Ties into `PLAN.md` §P2
(validation set) and §P5 (compare vs D3TaLES).

Two roles:
- **(A) Redox-annotated** — trusted oxidation/reduction potentials for validation,
  calibration (`E_exp = a·E_calc + b`), and as vetted candidates.
- **(B) Candidate pools** — large libraries with no redox annotation, filterable by
  redox motif (SMARTS/substructure search) to grow the screening set.

All URLs verified live (Aug 2026) via fetch/search. Publisher landing pages (RSC/ACS/
Elsevier) 403 automated fetch but resolve in-browser; the *data* lives on the repos below.

---

## Local status (downloaded) — as of 2026-08-20

Datasets live under `data/raw/validation/<name>/` and are **git-ignored** (bulk data is
never committed; re-clone/re-download with the commands below).

| Dataset | Local dir | State | Rows | Redox data |
|---|---|---|---|---|
| **OROP / ROP313** | `OROP/` | ✅ downloaded | 313 | **experimental** MeCN+DMF ox/red (V vs Fc) |
| **ReSolvedDB** | `ReSolvedDB/` | ✅ downloaded | 19,785 | computed reduction pot. (5 solvents incl. MeCN) |
| **D3TaLES** (bulk CC-BY) | `D3TaLES/d3tales_public.csv` | ✅ downloaded (35 MB) | 35,729 | computed ox/red (implicit MeCN, ε=35.688), λ, HOMO/LUMO, SA |
| **D3TaLES API** (code) | `D3TaLES_api/` | ✅ cloned | — | REST/processors/CV calculators (MIT) |

**Reproduce the downloads:**
```bash
cd data/raw/validation
git clone https://github.com/Liu-group/SI_data_redox_paper OROP     # experimental anchor
git clone https://github.com/grynova-ccc/ReSolvedDB                  # computed MeCN reductions
git clone https://github.com/D3TaLES/d3tales_api D3TaLES_api         # API code (no data)
# D3TaLES bulk dump: hosted on MDF/Globus (endpoint 82f1b5c6-...:/mdf_open/9bb8800f-.../1.1/).
# The HTTPS mirror serves known files but gives no directory listing, so list via authed
# globus-cli, then curl the file straight off the mirror (no Globus transfer needed):
globus ls "82f1b5c6-6e9b-11e5-ba47-22000b92c6ec:/mdf_open/9bb8800f-9cda-4957-ad12-60ad2a381177/1.1/"
mkdir -p D3TaLES && curl -L -o D3TaLES/d3tales_public.csv \
  "https://data.materialsdatafacility.org/mdf_open/9bb8800f-9cda-4957-ad12-60ad2a381177/1.1/d3tales_public.csv"
```

**D3TaLES utilization findings** (`scripts/analyze_d3tales.py` → `results/d3tales_*.csv`):
- Potentials are **absolute** (eV): `solv_oxidation_potential` −0.49..3.58, `solv_reduction_potential`
  3.42..11.00 → subtract a level-matched Fc absolute for V vs Fc. Ranking cancels the reference.
- Our-family coverage is **thin** (mostly ZINC drug-like neutrals): pyridinium 23, quinone 77,
  phenothiazine 23, **viologen 0**, metallocene 0 — but **nitroxide 990** (direct TEMPO-family probe).
- README **independently corroborates our viologen diagnosis**: it flags +2/−2 species as
  implicit-solvent-unreliable ("strong interactions with individual solvent molecules") — the
  same charge-dependent solvation failure we found for MV²⁺.
- **Best use:** (a) large *computed* cross-check for ranking/offset on overlapping chemistry;
  (b) the 990-nitroxide subset as a per-family validation probe; (c) λ (reorg.) + SA-score for
  multi-objective candidate ranking. **Not** a viologen/candidate source. Experimental gating
  stays with **OROP** (313 exp MeCN).

---

## Ranked table (best match first)

Rank = chemistry match × redox-data trust × ease of programmatic access.

| # | Database | Cat | Size | Redox data | Solvent | Access | License | Chemistry fit |
|---|----------|-----|------|-----------|---------|--------|---------|---------------|
| 1 | **D3TaLES** | A | ~35k–43k redox-active organics | ox/red potential, HOMO/LUMO, λ (reorg.), solubility, stability; **computed + experimental + literature** | **MeCN** (implicit, ε=35.688) | Python `d3tales_api` + REST (login); **CC-BY bulk dump on MDF/Globus, no login** | Data CC-BY 4.0, code MIT | **All 5 families** — purpose-built for non-aqueous RFB redoxmers |
| 2 | **ReSolvedDB** | A | 19,785 molecules / ~20k potentials | **reduction potential + EA**, computed | **MeCN** (SMD) + water/THF/DMSO/DMF | `git clone` GitHub (bulk) | CC-BY 4.0 | Strong — nitroxides, quinones, phenoxyl/phenazine radicals |
| 3 | **OROP / ROP313** | A | 313 (~183 organic) | **experimental** ox + red potential | **MeCN** + DMF | `git clone` GitHub (SI data) | repo (SI, permissive—verify) | Diverse organics; best **experimental MeCN anchor** for calibration |
| 4 | **MPcules** (Materials Project molecules) | A | ~170k molecules | ox/red free energy, IE, EA (**SHE-ref**) | vacuum/PCM/**SMD**; emphasis carbonate/Li — **MeCN unconfirmed** | `mp-api` / `pymatgen` (free API key) | CC-BY 4.0 | Good for N-heterocycles, quinones; open-shell via spin field |
| 5 | **Er/Aspuru-Guzik quinone screen (2015)** | A | ~1,710 quinone couples | computed reduction potential | **MeCN** (C-PCM) | RSC **ESI tables only** (no bulk repo) | RSC ESI | Quinones-in-MeCN, on-target but narrow + ESI-locked |
| 6 | **RedDB** | A | 31,618 (quinones + aza-aromatics) | redox potential (2H⁺/2e⁻), solubility, HOMO/LUMO, Fukui | **aqueous only** (PBF) | Harvard Dataverse (SQL/XLSX/CSV) | CC-BY 4.0 | Quinones + pyridinium-like aza-aromatics; **wrong solvent** |
| 7 | **OMol25** (Meta FAIR Chemistry) | B | ~83M systems / 100M+ SPs | **none** (energies+forces only; derive redox yourself) | gas-phase + **explicit** solvation incl. MeCN | HuggingFace (gated) + `fairchem`; Globus/MDF | Dataset CC-BY 4.0 | All families representable; all charge/spin states |
| 8 | **PubChem** | B | 110M+ compounds | none | n/a | PUG-REST / `pubchempy` (no key) | Public domain | Motif substructure search across all families |
| 9 | **ZINC22** | B | 230M+ purchasable | none | n/a | CartBlanche REST (no key) | free academic | Purchasable analogs by SMILES sim/substructure |
| 10 | **org-redox-dataset** (Kichev 2023) | A | quinones | computed redox (Li-reduction) | Li-electrode context (not MeCN) | `git clone` GitHub | GPL-3.0 | Quinones only; wrong measurement context |
| 11 | **ChEMBL** | B | 2.4M bioactives | none (bioactivity) | n/a | `chembl_webresource_client` | CC-BY-SA 3.0 | Low — bioactivity-focused, weak redox coverage |

---

## Per-database access recipes

### 1. D3TaLES — top pick (validation + candidates)
Digital Design to Transform Liquid-based Energy Storage — Risko group, Univ. of Kentucky
(Center for Applied Energy Research).

- Homepage: https://d3tales.as.uky.edu/ · Browse: https://d3tales.as.uky.edu/database/
- API docs: https://d3tales.github.io/d3tales_api/ · Code: https://github.com/D3TaLES/d3tales_api
  (MIT); DFT workflow https://github.com/D3TaLES/d3tales_fw ; schema https://github.com/D3TaLES/schema
- **Bulk CC-BY dump (NO login)** — MDF/Globus, DOI 10.18126/v5sj-6q93:
  https://acdc.alcf.anl.gov/mdf/detail/9bb8800f-9cda-4957-ad12-60ad2a381177-1.1/
- Paper: Duke, Risko et al., *Digital Discovery* **2**, 1152 (2023), DOI 10.1039/D3DD00081H
  (open: https://par.nsf.gov/servlets/purl/10508099)

**Programmatic query:**
```bash
pip install d3tales_api --use-deprecated=legacy-resolver   # pymatgen via conda
```
```python
from d3tales_api.D3database.restapi import D3talesData
# computed oxidation potential across the DB -> pandas DataFrame
df = D3talesData().get_prop_data('mol_characterization.oxidation_potential.value')
# also: reduction_potential, solvation_energy, reorganization_energy, HOMO/LUMO ...
```
- **Auth:** REST API needs a login (cookie + CSRF POST; creds via `UPLOAD_USER`/`UPLOAD_PASS`
  env vars). Account signup is not documented online — register on the site or contact the
  Risko group. **Public web browse and the MDF bulk dump need no account** → easiest entry.
- **Size:** ~35,729 public in the dump / ~38,434 advertised for ML / >43,000 total.
- **Properties:** >90 per molecule (eV). ox/red potential, HOMO/LUMO, **reorganization
  energy**, solubility, stability. `solv_*` props are in **implicit acetonitrile (ε=35.688)**.
  Check the Digital Discovery SI for the exact reference-electrode convention.
- **Theory / provenance:** DFT **(IP-tuned) LC-ωHPBE/Def2SVP** (Gaussian16, Fireworks HT);
  plus **experimental** uploads (CV/robotic) and **literature NLP** (ChemDataExtractor2).
- **Fit:** direct — built for non-aqueous MeCN organic RFB; spans all five families.
  Use the API/SMILES filter to pull per-family subsets for §P2/§P5 number-to-number checks.

### 2. ReSolvedDB — best downloadable computed MeCN set
Multi-solvent GNN reduction-potential dataset.

- Dataset: https://github.com/grynova-ccc/ReSolvedDB · Models: https://github.com/grynova-ccc/ReSolved
- Paper: *J. Chem. Inf. Model.* 2025, DOI 10.1021/acs.jcim.5c01450
- **Get it:** `git clone https://github.com/grynova-ccc/ReSolvedDB`
- **Size:** 19,785 molecules; ~20k reduction potentials over 5 solvents.
- **Properties:** reduction potential + electron affinity in **acetonitrile**, water, THF,
  DMSO, DMF. SMILES + DFT geometries included. Computed (not experimental).
- **Theory:** **M06-2X/def2-TZVPD, SMD** implicit solvent.
- **License:** CC-BY 4.0.
- **Fit:** strong — nitroxides, phenoxyl radicals, quinones, phenazines; closed/open-shell
  radicals. Use as a computed MeCN cross-check and candidate expansion.

### 3. OROP / ROP313 — experimental MeCN anchor (calibration ground truth)
- Repo: https://github.com/Liu-group/SI_data_redox_paper
- Papers: *J. Chem. Theory Comput.* 2022, DOI 10.1021/acs.jctc.1c01040; review 2023
  DOI 10.1021/acs.jctc.3c00355 (PMC10414033)
- **Get it:** `git clone https://github.com/Liu-group/SI_data_redox_paper`
- **Size:** ROP313 = 313 organic/organometallic; OROP ≈ 183–193 purely organic.
- **Properties:** **experimentally measured** reduction AND oxidation potentials in
  **acetonitrile** + DMF; repo also ships GFN2-xTB and PWPB95-D4 computed values.
- **Fit:** the curated experimental MeCN set used for ML benchmarking — ideal for the
  §V `E_exp = a·E_calc + b` calibration + MAE. Elements H,C,N,O,S,F,Cl,Br,I,Na,B.
- **License:** SI-data repo, not explicitly stated — verify before redistribution.

### 4. MPcules — Materials Project molecules (Electrolyte Genome heritage)
Persson/Blau groups, LBNL.

- Docs: https://docs.materialsproject.org/methodology/molecules-methodology/redox-and-electrochemical-properties
- Paper: Spotte-Smith et al., *Digital Discovery* **2**, 1862 (2023), DOI 10.1039/D3DD00153A
- **API key:** free account at https://next-gen.materialsproject.org → key at `/api`.
```bash
pip install mp-api pymatgen
```
```python
from mp_api.client import MPRester
with MPRester("YOUR_API_KEY") as mpr:
    docs = mpr.molecules.summary.search(
        elements=["C","N"], nelements=(2,30),
        fields=["molecule_id","formula_alphabetical","charge","spin_multiplicity",
                "ionization_energy","electron_affinity",
                "reduction_free_energy","oxidation_free_energy"])
```
- **Size:** ~170k DFT molecules. **Redox:** IE, EA, adiabatic ox/red free energies,
  **SHE-referenced** ox/red potentials (from `MoleculeThermoDoc`s).
- **Theory:** ωB97X-D/ωB97X-V/ωB97M-V; def2-SVPD/TZVPPD/QZVPPD; PCM or SMD (QChem).
- **Solvent caveat:** heavy carbonate/Li-ion (LIBE/MADEIRA) emphasis; **MeCN not confirmed** —
  no top-level `solvent=` filter, inspect the solvent tag on returned redox docs.
- **License:** CC-BY 4.0. Good IE/EA cross-check; open-shell radicals via `spin_multiplicity`.

### 5. Er / Aspuru-Guzik all-quinone RFB screen (2015)
- Paper: Er, Suh, Marshak, Aspuru-Guzik, *Chem. Sci.* **6**, 885 (2015), DOI 10.1039/C4SC03030C
  (free full text: https://www.osti.gov/servlets/purl/1624913)
- ~1,710 quinone/hydroquinone couples; computed reduction potentials in **acetonitrile**
  (B3LYP-class, **C-PCM** MeCN). **Access limitation:** data is in the RSC **ESI tables**
  (PDF/spreadsheet), no Dataverse/figshare mirror → parse the ESI. On-target but narrow.
- Note: the broader **Harvard Clean Energy Project (CEPDB, molecularspace.org)** is an
  organic-**photovoltaic** database (HOMO/LUMO/PCE), **not** redox-in-MeCN — not a match.

### 6. RedDB — aqueous RFB database (complementary)
- Paper: Sorkun et al., *Sci. Data* **9**, 718 (2022), DOI 10.1038/s41597-022-01832-2
- Download: Harvard Dataverse DOI 10.7910/DVN/F3QFSQ
  (https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/F3QFSQ);
  code https://github.com/ergroup/RedDB . Formats: SQL, XLSX, CSV.
- 31,618 molecules (quinones + aza-aromatics), redox potential + ML solubility + HOMO/LUMO +
  Fukui. **Theory:** PBE/LACVP**, **aqueous (PBF)** only → **wrong solvent**; quinone/
  pyridinium chemistry usable for transfer-learning, not direct MeCN comparison. CC-BY 4.0.

### 7. OMol25 — Meta FAIR Chemistry (candidate pool / MLIP data)
- HuggingFace (gated): https://huggingface.co/facebook/OMol25 · Docs: https://fair-chem.github.io/omol25/
- Paper: Levine et al., arXiv:2505.08762 (2025). **This is the OMol source behind UMA** — the
  same theory the project's UMA pre-optimizer was trained on.
- **Access:** accept HF terms (legal name/DOB/org) → `huggingface-cli login`; load via
  `fairchem`:
```bash
pip install fairchem-core
```
```python
from fairchem.core.datasets import AseDBDataset
ds = AseDBDataset({"src": "/path/to/omol25/train"})   # *.aselmdb
atoms = ds.get_atoms(0); atoms.info["charge"], atoms.info["spin_multiplicity"]
```
- **Size:** ~83M systems / 100M+ single-points, ≤350 atoms, 83 elements.
- **Properties:** total energy + forces **only** — **no redox potentials** (derive IE/EA/ΔG
  from same-connectivity charge/spin pairs yourself). Domains include an **electrolytes**
  subset with **explicit** solvation incl. **MeCN**; all charge/spin states present.
- **Theory:** **ωB97M-V/def2-TZVPD** (ORCA6). Dataset **CC-BY 4.0** (model weights: separate
  FAIR Chemistry License). Best as MLIP training data or a charged/open-shell candidate pool.

### 8–9. PubChem & ZINC22 — motif substructure candidate pools (Category B)
Use validated SMARTS/SMILES below to pull motif-matched candidates.

**PubChem** (`pip install pubchempy`, no key; ≤5 req/s):
```python
import pubchempy as pcp
hits = pcp.get_compounds('c1ccc2c(c1)Nc1ccccc1S2', 'smiles',   # phenothiazine core
                         searchtype='substructure', MaxRecords=500)
```
**ZINC22** (CartBlanche REST, no key) — purchasable analogs by SMILES:
```bash
curl "https://cartblanche22.docking.org/smiles.txt:smiles=O=C1C=CC(=O)C=C1&dist=4&output_fields=zinc_id,smiles,catalogs"
```
**ChEMBL** (`pip install chembl_webresource_client`) — `new_client.substructure.filter(smiles=...)`;
lowest priority (bioactivity, not redox).

**Validated redox-motif patterns** (RDKit-checked) for substructure/SMARTS filtering:

| Motif | Pattern |
|-------|---------|
| pyridinium | `[n+]1ccccc1` |
| viologen (4,4′-bipyridinium) | `[n+]1ccc(-c2cc[n+](cc2))cc1` |
| quinone (para) | `O=C1C=CC(=O)C=C1` |
| anthraquinone | `O=C1c2ccccc2C(=O)c2ccccc21` |
| phenothiazine | `c1ccc2c(c1)Nc1ccccc1S2` |
| nitroxide / aminoxyl | `[#6][N]([#6])[O]` |
| TEMPO core | `CC1(C)CCCC(C)(C)N1[O]` |

---

## Recommendation — integrate these first

1. **D3TaLES** *(direct match; validation + candidates)*.
   First step — grab the CC-BY bulk dump (no account) and index by SMILES:
   download the MDF/Globus package (DOI 10.18126/v5sj-6q93,
   https://acdc.alcf.anl.gov/mdf/detail/9bb8800f-9cda-4957-ad12-60ad2a381177-1.1/), then
   filter molecules by the motif SMARTS above and pull `oxidation_potential`/
   `reduction_potential` + `solv_*` (MeCN) columns. In parallel `pip install d3tales_api`
   and try `D3talesData().get_prop_data('mol_characterization.oxidation_potential.value')`;
   if the REST API demands a login, stay on the dump. Feeds §P5 (compare our E° to theirs).

2. **OROP/ROP313** *(experimental MeCN ground truth for calibration)*.
   First step — `git clone https://github.com/Liu-group/SI_data_redox_paper`, extract the
   organic MeCN oxidation/reduction potentials, and use them as the `E_exp` side of the
   §V `E_exp = a·E_calc + b` fit / MAE. Small, experimental, directly on-solvent.

3. **ReSolvedDB** *(computed MeCN cross-check + candidate expansion)*.
   First step — `git clone https://github.com/grynova-ccc/ReSolvedDB`, keep the acetonitrile
   reduction potentials + SMILES, and use the nitroxide/quinone/radical entries both to
   sanity-check our M06-2X-class numbers and to widen the candidate list.

**Also worth wiring in early:** MPcules via `mp-api` (free key) for an independent IE/EA/
SHE-referenced redox cross-check — verify the solvent tag per doc since MeCN coverage is
unconfirmed. Treat **OMol25** as the charged/open-shell candidate pool + MLIP data (it is the
UMA training source, so it aligns with the project's pre-optimizer), and **RedDB** as
aqueous quinone/pyridinium transfer-learning data only.
