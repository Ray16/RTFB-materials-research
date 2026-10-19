"""External optimizer scores retain their actual ordinal meaning."""
import csv
import hashlib
import json

import pytest

from redox.screening.source_audit import audit_sources, optimizer_grade, reaction_quality


def test_zero_is_the_best_source_grade_and_four_is_nonoptimal():
    assert reaction_quality("0", "0.0") == "both_best_optimizer_grade"
    assert reaction_quality("0", "4") == "nonoptimal_geometry"
    assert reaction_quality("4", "0") == "nonoptimal_geometry"
    assert reaction_quality("1", "3") == "intermediate_grade_review"


@pytest.mark.parametrize("value", [None, "", True, False, "True", "nan", "inf", "-1", "5", "2.5"])
def test_malformed_or_boolean_source_grade_cannot_look_valid(value):
    with pytest.raises(ValueError):
        optimizer_grade(value)
    assert reaction_quality(value, 0) == "invalid_or_missing_grade"


def test_audit_retains_identifiers_and_raw_flags_without_admitting_labels(tmp_path):
    source = tmp_path / "raw"; source.mkdir()
    reaction = source / "RedDBv2_reaction.tab"
    reaction.write_text("reaction_id\treactant_convergence\tproduct_convergence\n"
                        "A\t0\t0\nB\t0\t4\nC\t2\t3\nD\t\t0\n")
    resolved = source / "ReSolvedData.csv"
    resolved.write_text("same_sml\nyes\nsame\nno\nidk\n")
    before = {p: p.read_bytes() for p in (reaction, resolved)}
    output = tmp_path / "audit"
    report = audit_sources(source, output)
    assert report["reddb_quality_counts"] == {"both_best_optimizer_grade": 1,
        "nonoptimal_geometry": 1, "intermediate_grade_review": 1, "invalid_or_missing_grade": 1}
    assert report["resolved_same_sml_counts"] == {"yes": 1, "same": 1, "no": 1, "idk": 1}
    assert report["source_labels_admitted"] == report["quantum_jobs_submitted"] == 0
    rows = list(csv.DictReader((output / "reddb_reaction_quality.csv").open()))
    assert [r["reaction_id"] for r in rows] == ["A", "B", "C", "D"]
    assert [int(r["source_row"]) for r in rows] == [2, 3, 4, 5]
    assert all(r["project_label_admission"].startswith("not_admitted") for r in rows)
    for path, content in before.items():
        assert path.read_bytes() == content
        assert report["source_sha256"][str(path.resolve())] == hashlib.sha256(content).hexdigest()
    assert json.loads((output / "summary.json").read_text()) == report
    with pytest.raises(FileExistsError):
        audit_sources(source, output)
