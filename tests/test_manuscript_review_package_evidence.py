"""The manuscript must distinguish cached replay, public release and fresh QC."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("release_claims", ROOT / "manuscript/check_claims.py")
CLAIMS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLAIMS)
EPOCH = ROOT / "manuscript/evidence/2026-10-23_TOC_release_revision"


def load(name):
    return json.loads((EPOCH / name).read_text())


def test_frozen_release_coverage_reconstructs_without_promoting_scope():
    result = CLAIMS.audit_reviewer_package(load)
    assert result == dict(raw_states=243, electronic_only_records=32, archive_files=2498,
                          Hessian_states=242, negative_checks=7,
                          public_deposition=False, full_GPU_restoration=False)


@pytest.mark.parametrize("filename,change,error", [
    ("review_package_manifest.json", lambda x: x.update(raw_states=244), "inventory counts"),
    ("review_package_data.json", lambda x: x["states"][0].update(partition_kind="monatomic"), "Hessian coverage"),
    ("review_package.json", lambda x: x.update(public_archive_deposited=True), "deposition scope"),
    ("review_package_replay.json", lambda x: x.update(full_GPU_environment_restored=True), "deposition scope"),
    ("review_package_hessians.json", lambda x: x.update(original_physical_QC_repeated=True), "fresh physical"),
    ("review_package_arithmetic.json", lambda x: x.update(benchmark_endpoints=38), "primary replay coverage"),
    ("review_package_replay.json", lambda x: x["negative_checks"][0].update(rejected=False), "negative checks"),
])
def test_inflated_or_relabelled_release_claims_are_rejected(filename, change, error):
    def altered(name):
        value = load(name)
        if name == filename:
            change(value)
        return value
    with pytest.raises(ValueError, match=error):
        CLAIMS.audit_reviewer_package(altered)
