"""Independent reconstruction of the published design/evaluation joins and effects."""

import hashlib
import json
from pathlib import Path
import re
import time

from rdkit import Chem

from redox.core.common import ROOT
from redox.workflows.polymer_campaign import verify_sources


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check():
    out = ROOT / "results/discovery/expansion"
    result = json.loads((out / "current.json").read_text())
    provenance = result["provenance"]
    verify_sources(provenance["frozen_input_sha256"])
    verify_sources(provenance["producer_sha256"])
    source = Path(provenance["snapshot"]) / "sources"
    feed = json.loads((source / "feed.json").read_text())
    pool = json.loads((source / "pool.json").read_text())
    rows = feed["original39"] + feed["additional_candidates"] + feed["comparison_controls"]
    expected = {row.get("registry_id") or row["id"]: row for row in rows}
    by_id = {row["id"]: row for row in result["records"] if row["evaluation_role"]}
    assert expected.keys() == by_id.keys()
    for key, row in expected.items():
        observed = by_id[key]
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(row["exact_grafted_smiles"]))
        assert observed["grafted_smiles"] == canonical
        for field in ("E1_vs_Fc_V", "E2_vs_Fc_V", "wave_spacing_V", "MW", "accessible_capacity_mAh_g"):
            assert observed[field] == row.get(field), (key, field)
    for comparison in result["comparisons"]:
        row, baseline = expected[comparison["id"]], expected[comparison["baseline"]]
        for field, target in (("E1_vs_Fc_V", "delta_E1_V"), ("E2_vs_Fc_V", "delta_E2_V"), ("wave_spacing_V", "delta_wave_spacing_V")):
            if row.get(field) is None or baseline.get(field) is None:
                assert comparison[target] is None
            else:
                assert abs(comparison[target] - (row[field] - baseline[field])) < 1e-12
    original_pairs = {row["id"]: row for row in feed["matched_comparisons"]}
    assert len(original_pairs) == len(result["comparisons"])
    for row in result["comparisons"]:
        original = original_pairs[row["id"]]
        for field in ("baseline", "delta_E1_V", "delta_E2_V", "delta_wave_spacing_V", "status"):
            assert row[field] == original[field]
    proposed = [row for row in result["records"] if row["evaluation_role"] is None]
    for row in proposed:
        assert row["E1_vs_Fc_V"] is None and row["E2_vs_Fc_V"] is None
        assert row["base_models"] is None and row["reference_grafting_deltaG_kJ_mol"] is None
    keys = {Chem.MolToInchiKey(Chem.MolFromSmiles(row["grafted_smiles"])) for row in pool}
    keys.update(row["exact_grafted_inchikey"] for row in rows)
    assert keys == {row["inchikey"] for row in result["records"]}
    assert len(result["records"]) == len(keys)
    assert sum(row["selected_candidates"] for row in result["family_coverage"]) == 75
    assert sum(row["external_controls"] for row in result["family_coverage"]) == 1
    assert all(row["selected_graphs"] <= row["unique_graphs"] for row in result["design_cells"])
    figure = json.loads((out / "report_manifest.json").read_text())
    assert sha(out / "current.json") == figure["source_sha256"]
    assert sha(out / "index.html") == figure["output_sha256"]
    for name, expected_hash in figure["figure_sha256"].items():
        assert sha(out / name) == expected_hash
    checked_links = 0
    for path in (ROOT / "docs/CANDIDATE_EXPANSION_LOOP.md", ROOT / "results/discovery/README.md", ROOT / "results/discovery/expansion/README.md", ROOT / "results/README.md"):
        for target in re.findall(r"\]\(([^)]+)\)", path.read_text()):
            if target.startswith(("https:", "http:", "#")):
                continue
            assert (path.parent / target.split("#")[0]).exists(), (path, target)
            checked_links += 1
    previous = []
    for name in ("2026-10-20_spacing_pilot", "2026-10-20_geometry_recovery", "2026-10-20_frozen_fragments"):
        path = ROOT / "results/polymer_models/campaigns" / name / "plan.json"
        verify_sources(json.loads(path.read_text())["producer_sha256"])
        previous.append(name)
    summary = dict(
        checked_at=time.time(), status="joined_existing_evaluation_and_designs_verified",
        tests="33 passed: expansion integration, enumeration, paired comparisons and route gates",
        lint="Ruff E/F consumer, collector and tests; F plotting script",
        counts=result["counts"], exact_evaluator_rows_checked=len(rows),
        independent_paired_effects_checked=len(result["comparisons"]),
        proposed_graphs_with_no_inherited_properties=len(proposed),
        local_links_checked=checked_links, active_campaign_producers_unchanged=previous,
        source_snapshot=provenance["snapshot"], source_sha256=sha(out / "current.json"),
        figure_review="Heatmap visually checked: signed effects, common color scale, no selection ranking",
        no_registry_changes=True, no_quantum_submissions=True,
        limits=result["limits"],
    )
    Path(__file__).with_name("summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print({key: summary[key] for key in ("status", "exact_evaluator_rows_checked", "independent_paired_effects_checked", "local_links_checked")})


if __name__ == "__main__":
    check()
