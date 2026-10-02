"""Fill data/raw/candidates/Candidates.xlsx with the values computed by this pipeline.

The sheet's "SMILES string (neutral monosubstituted)" column is the STANDALONE (pre-grafting)
molecule, so each row maps to a standalone library id; the resin-grafted analogue that the
project actually targets is reported alongside in `*_grafted` columns.

Hole vs electron reorganization is keyed on n_e from the manifest, NOT on net charge:
`ndi_ammonium` carries a permanent +1 trimethylammonium, so its first reduction runs
+1 -> 0 and would be misread as an oxidation by charge alone.

  lambda_electron = first reduction  (n_e:  0 -> -1)
  lambda_hole     = first oxidation  (n_e: +1 <-  0)
  lambda_e2       = second reduction (n_e: -1 -> -2)

Original sheet columns are left untouched (they hold the user's literature / D3TaLES values);
everything computed here lands in new, explicitly-named columns.

QC: a lambda that redox.properties.reorg flagged (conformer_jump / anion_unbound / negative_*) is NOT
written to the clean `lambda_*_eV` column — it goes to `lambda_*_eV_unfiltered` alongside the
flag and the heavy-atom RMSD, so a contaminated value cannot be read as a validated one.

  python scripts/pipeline/fill_candidates_xlsx.py
"""
from __future__ import annotations
import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
XLSX_IN = ROOT / "data" / "raw" / "candidates" / "Candidates.xlsx"
XLSX_OUT = RESULTS / "Candidates_computed.xlsx"
CSV_OUT = RESULTS / "starting_candidates_summary.csv"

# sheet row order -> (standalone id matching the sheet's SMILES, resin-grafted analogue)
ROW_MAP = [
    ("pmdi_sa",          "pmdi"),           # sheet: N-2-pentyl PMDI; standalone modelled as N,N'-dimethyl
    ("methyl_viologen",  "viologen"),
    ("ethylviologen_sa", "ethylviologen"),
    ("ndi_ammonium_sa",  "ndi_ammonium"),
    ("mophquinone_sa",   "mophquinone"),
    ("dmophquinone_sa",  "dmophquinone"),
]


def _read(name):
    with (RESULTS / name).open() as f:
        return list(csv.DictReader(f))


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def couple_labels(manifest_rows):
    """id -> {couple_string: role} where role is lambda_electron / lambda_hole / lambda_e2."""
    by_id = {}
    for r in manifest_rows:
        by_id.setdefault(r["id"], {})[r["state"]] = int(r["n_e"])
    out = {}
    for gid, states in by_id.items():
        ne2state = {v: k for k, v in states.items()}
        roles = {}
        if 0 in ne2state and -1 in ne2state:
            roles[f"{ne2state[0]}->{ne2state[-1]}"] = "lambda_electron"
        if 1 in ne2state and 0 in ne2state:
            roles[f"{ne2state[1]}->{ne2state[0]}"] = "lambda_hole"
        if -1 in ne2state and -2 in ne2state:
            roles[f"{ne2state[-1]}->{ne2state[-2]}"] = "lambda_e2"
        out[gid] = roles
    return out


def collect(gid, reorg, pots, cap, roles):
    """All computed quantities for one library id."""
    d = {}
    for r in reorg:
        if r["id"] != gid:
            continue
        role = roles.get(gid, {}).get(r["couple"])
        if not role:
            continue
        # A lambda carrying a QC flag from redox.properties.reorg is not a usable inner-sphere value.
        # Report it in a separate *_unfiltered column with the flag, and leave the clean
        # column empty, so a flagged number can never be mistaken for a validated one.
        qc = (r.get("flag") or "").strip()
        val = _f(r["lambda_i_eV"])
        d[f"{role}_eV_unfiltered"] = val
        d[f"{role}_qc_flag"] = qc
        d[f"{role}_rmsd_A"] = _f(r.get("rmsd_A"))
        d[f"{role}_eV"] = None if qc else val
    for r in pots:
        if r["id"] != gid:
            continue
        role = roles.get(gid, {}).get(r["event"])
        if role == "lambda_electron":
            d["E_first_reduction_V_vs_Fc"] = _f(r["E_vs_Fc_V"])
        elif role == "lambda_hole":
            d["E_first_oxidation_V_vs_Fc"] = _f(r["E_vs_Fc_V"])
        elif role == "lambda_e2":
            d["E_second_reduction_V_vs_Fc"] = _f(r["E_vs_Fc_V"])
    for r in cap:
        if r["id"] == gid:
            d["SA_score"] = _f(r["SA_score"])
            d["MW"] = _f(r["MW"])
            d["specific_capacity_mAh_g"] = _f(r["specific_capacity_mAh_g"])
            d["n_electrons"] = _f(r["n_electrons"])
            d["dGsolv_neutral_eV"] = _f(r["dGsolv_neutral_eV"])
    return d


FIELDS = ["SA_score",
          "lambda_electron_eV", "lambda_e2_eV", "lambda_hole_eV",
          "lambda_electron_eV_unfiltered", "lambda_e2_eV_unfiltered",
          "lambda_hole_eV_unfiltered",
          "lambda_electron_qc_flag", "lambda_e2_qc_flag", "lambda_hole_qc_flag",
          "lambda_electron_rmsd_A", "lambda_e2_rmsd_A", "lambda_hole_rmsd_A",
          "E_first_reduction_V_vs_Fc", "E_second_reduction_V_vs_Fc",
          "E_first_oxidation_V_vs_Fc", "n_electrons", "MW",
          "specific_capacity_mAh_g", "dGsolv_neutral_eV"]


def main():
    manifest = list(csv.DictReader((ROOT / "library" / "manifest.csv").open()))
    roles = couple_labels(manifest)
    reorg, pots, cap = _read("reorganization.csv"), _read("redox_potentials.csv"), \
        _read("capacity_and_proxies.csv")

    df = pd.read_excel(XLSX_IN)
    if len(df) != len(ROW_MAP):
        raise SystemExit(f"sheet has {len(df)} rows, ROW_MAP has {len(ROW_MAP)} — remap first")

    flat = []
    for i, (sa_id, graft_id) in enumerate(ROW_MAP):
        s = collect(sa_id, reorg, pots, cap, roles)
        g = collect(graft_id, reorg, pots, cap, roles)
        df.loc[i, "computed_id_standalone"] = sa_id
        df.loc[i, "computed_id_grafted"] = graft_id
        for k in FIELDS:
            if s.get(k) is not None:
                df.loc[i, f"{k}_computed"] = s[k]
            if g.get(k) is not None:
                df.loc[i, f"{k}_grafted"] = g[k]
        flat.append(dict(sheet_row=i, candidate=str(df.iloc[i, 0]).split("\n")[0],
                         id_standalone=sa_id, id_grafted=graft_id,
                         **{f"{k}_standalone": s.get(k) for k in FIELDS},
                         **{f"{k}_grafted": g.get(k) for k in FIELDS}))

    df.to_excel(XLSX_OUT, index=False)
    pd.DataFrame(flat).to_csv(CSV_OUT, index=False)
    print(f"[xlsx] {XLSX_OUT.relative_to(ROOT)}  ({len(df)} rows, "
          f"{len(df.columns)} cols)")
    print(f"[csv ] {CSV_OUT.relative_to(ROOT)}")
    show = ["candidate", "id_standalone", "SA_score_standalone",
            "lambda_electron_eV_standalone", "lambda_electron_qc_flag_standalone",
            "lambda_electron_eV_grafted", "lambda_electron_qc_flag_grafted"]
    print()
    print(pd.DataFrame(flat)[show].to_string(index=False))


if __name__ == "__main__":
    main()
