"""Audit downloaded source quality without admitting labels or registering molecules."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path

from redox.core.io import atomic_text_writer, atomic_write_json

REDDB_DOCUMENTATION = "https://www.nature.com/articles/s41597-022-01832-2"


def optimizer_grade(value):
    """RedDB uses ordinal Jaguar grades 0–4, not Boolean convergence."""
    if isinstance(value, bool):
        raise ValueError("Boolean convergence is not a RedDB optimizer grade")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("Missing or invalid optimizer grade") from exc
    if not math.isfinite(number) or number not in range(5):
        raise ValueError("Optimizer grade must be an integer from zero through four")
    return int(number)


def reaction_quality(reactant, product):
    try:
        grades = optimizer_grade(reactant), optimizer_grade(product)
    except ValueError:
        return "invalid_or_missing_grade"
    if 4 in grades:
        return "nonoptimal_geometry"
    if grades == (0, 0):
        return "both_best_optimizer_grade"
    return "intermediate_grade_review"


def audit_sources(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    reddb, resolved = source / "RedDBv2_reaction.tab", source / "ReSolvedData.csv"
    fingerprints = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (reddb, resolved)}
    with reddb.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows:
        raise ValueError("Source reaction table is empty")
    quality = [dict(source_row=i, reaction_id=r["reaction_id"],
                    reactant_convergence=r["reactant_convergence"],
                    product_convergence=r["product_convergence"],
                    quality=reaction_quality(r["reactant_convergence"], r["product_convergence"]),
                    project_label_admission="not_admitted; aqueous_PCET_source_model_only")
               for i, r in enumerate(rows, 2)]
    with resolved.open(newline="") as stream:
        identity_flags = Counter(r.get("same_sml", "") for r in csv.DictReader(stream))
    output.mkdir(parents=True, exist_ok=False)
    with atomic_text_writer(output / "reddb_reaction_quality.csv", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(quality[0]))
        writer.writeheader(); writer.writerows(quality)
    report = dict(source_sha256=fingerprints, reddb_rows=len(rows),
                  reddb_quality_counts=dict(Counter(r["quality"] for r in quality)),
                  reddb_grade_documentation=REDDB_DOCUMENTATION,
                  reddb_policy="Grade 0 is best; 4 denotes a nonoptimal structure. Grades 1–3 retain review status; optimizer grade alone does not establish accuracy.",
                  resolved_same_sml_counts=dict(identity_flags),
                  resolved_policy="Raw flags preserved; yes/same/no/idk semantics not inferred from spelling. No source labels admitted.",
                  source_labels_admitted=0, quantum_jobs_submitted=0)
    changed = [str(p) for p in (reddb, resolved)
               if hashlib.sha256(p.read_bytes()).hexdigest() != fingerprints[str(p)]]
    if changed:
        raise RuntimeError(f"Source changed during quality audit: {changed}")
    atomic_write_json(output / "summary.json", report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="preserved directory containing both downloaded source tables")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(audit_sources(args.source, args.out), indent=2))


if __name__ == "__main__":
    main()
