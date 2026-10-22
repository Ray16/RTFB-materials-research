"""Frozen experimental panels and descriptive, dependency-aware error reports.

No calculation is admitted here: the independent matched-parent reviewer owns
state/reference/spin QC. This consumer verifies that review and its frozen inputs.
Historical measurements and model development are explicitly retrospective.
"""

from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def key(row):
    return tuple(str(row[name]) for name in ("id", "event", "doi", "source", "where"))


def number(value):
    if value in (None, ""):
        return None
    if isinstance(value, bool):
        raise ValueError("Boolean is not a benchmark value")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Nonfinite benchmark value")
    return result


def metrics(errors):
    errors = [number(e) for e in errors]
    if any(e is None for e in errors):
        raise ValueError("Missing residual cannot be scored")
    n = len(errors)
    return dict(
        n=n,
        MAE_V=sum(abs(e) for e in errors) / n if n else None,
        RMSE_V=math.sqrt(sum(e * e for e in errors) / n) if n else None,
        bias_V=sum(errors) / n if n else None,
        maximum_absolute_error_V=max(map(abs, errors)) if n else None,
    )


def family(description):
    return next(
        (
            name
            for name in (
                "quinone",
                "aromatic imide",
                "viologen",
                "nitroxide",
                "phenothiazine",
            )
            if description.startswith(name)
        ),
        description,
    )


def freeze(entries, sources):
    """Freeze strict endpoints before completion; never claim a holdout."""
    catalog, strict = [], []
    for molecule in entries:
        waves = Counter()
        for event in molecule["events"]:
            if event["kind"] != "E":
                continue
            source_key = (event["doi"], event["source"])
            waves[source_key] += 1
            if waves[source_key] > 2:
                raise ValueError("Expected at most two ordered waves per source")
            row = dict(
                id=molecule["id"],
                name=molecule["name"],
                parent_smiles=molecule["smiles"],
                family=family(molecule["bench_family"]),
                source_family=molecule["bench_family"],
                event=event["event"],
                wave=waves[source_key],
                exp_V=number(event["exp_V_vs_Fc"]),
                tier=event["tier"],
                source=event["source"],
                doi=event["doi"],
                where=event["where"],
                source_limitation=event["tier_reason"],
                excluded_molecule=molecule.get("excluded"),
            )
            if row["exp_V"] is None:
                raise ValueError("Experimental target absent")
            catalog.append(row)
            if row["tier"] == "A" and not row["excluded_molecule"]:
                if not row["doi"]:
                    raise ValueError("Strict experimental source absent")
                strict.append(row)
    if not strict or len({key(r) for r in catalog}) != len(catalog):
        raise ValueError("Empty strict panel or duplicate source endpoint")
    payload = dict(
        schema=1,
        scope=(
            "Exact free parents in sourced MeCN conditions; not exact-graft validation"
        ),
        study_design="retrospective; all currently sourced strict endpoints included",
        held_out_from_model_development=False,
        fitted_correction=False,
        endpoint_selection=(
            "Tier A and no molecule-level exclusion; "
            "independent computed QC still required"
        ),
        primary_analysis=(
            "Unrounded potential residuals; wave- and family-specific "
            "errors; all outliers retained"
        ),
        uncertainty_analysis=(
            "Molecule and original-study cluster bootstrap; "
            "descriptive of this panel only"
        ),
        minimum_clusters_for_interval=5,
        bootstrap_replicates=2000,
        bootstrap_seed=20261021,
        interval_percentiles=[2.5, 97.5],
        source_sha256={str(Path(p).resolve()): digest(p) for p in sources},
        targets=sorted(strict, key=key),
        source_catalog=sorted(catalog, key=key),
    )
    payload["panel_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()
    ).hexdigest()
    return payload


def validate_panel(panel):
    expected = panel["panel_sha256"]
    payload = {k: v for k, v in panel.items() if k != "panel_sha256"}
    if (
        hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
        != expected
    ):
        raise ValueError("Frozen panel content changed")
    verify_hashes(panel["source_sha256"])


def verify_hashes(inventory):
    if not inventory:
        raise ValueError("Evidence inventory is empty")
    for name, expected in inventory.items():
        if digest(name) != expected:
            raise ValueError(f"Changed benchmark evidence: {name}")


def cluster_intervals(rows, field, *, replicates=2000, seed=20261021, minimum=5):
    """Resample clusters, retaining every endpoint and multiplicity within a draw."""
    import numpy as np

    clusters = defaultdict(list)
    for row in rows:
        clusters[row[field]].append(row["error_V"])
    count = len(clusters)
    result = dict(cluster=field, independent_cluster_count=count, intervals=None)
    if count < minimum:
        result["status"] = "insufficient_clusters_for_interval"
        return result
    values = [clusters[name] for name in sorted(clusters)]
    rng = np.random.default_rng(seed)
    draws = []
    for indices in rng.integers(0, count, size=(replicates, count)):
        # Repeated sampled clusters must remain repeated, not be deduplicated.
        draws.append(metrics([error for i in indices for error in values[i]]))
    result.update(
        status="panel_resampling_only",
        replicates=replicates,
        seed=seed,
        intervals={
            name: np.percentile([d[name] for d in draws], [2.5, 97.5]).tolist()
            for name in ("MAE_V", "RMSE_V", "bias_V")
        },
    )
    return result


def summarize(panel, comparisons):
    targets = {key(r): r for r in panel["targets"]}
    if len(targets) != len(panel["targets"]):
        raise ValueError("Repeated target")
    by_key = {}
    for row in comparisons:
        if row["kind"] != "E":
            continue
        identity = key(row)
        if identity in by_key:
            raise ValueError("Repeated comparison endpoint")
        by_key[identity] = row
        if row["tier"] == "A" and identity not in targets:
            raise ValueError("Strict comparison outside the frozen panel")
    scored, pending = [], []
    for identity, target in targets.items():
        row = by_key.get(identity)
        if row is None:
            pending.append(dict(target, status="missing comparison"))
            continue
        if number(row["exp_V"]) != target["exp_V"] or row["tier"] != "A":
            raise ValueError("Experimental identity/target/tier changed")
        published = number(row["calc_V"])
        if published is None:
            pending.append(dict(target, status=row["status"]))
            continue
        value = number(row.get("calc_unrounded_V"))
        if value is None or row["status"] != "checked_native_parent":
            raise ValueError(
                "Calculated value lacks independent admission/unrounded evidence"
            )
        if abs(value - published) > 0.00050001:
            raise ValueError("Published and unrounded values disagree")
        serialized_error = number(row["error_V"])
        if (
            serialized_error is None
            or abs(published - target["exp_V"] - serialized_error) > 1e-9
        ):
            raise ValueError("Serialized benchmark residual disagrees")
        scored.append(
            dict(
                target,
                calc_V=value,
                serialized_calc_V=published,
                error_V=value - target["exp_V"],
            )
        )

    def group(field):
        return {
            str(value): dict(
                metrics([r["error_V"] for r in scored if r[field] == value]),
                molecules=len({r["id"] for r in scored if r[field] == value}),
                original_studies=len({r["doi"] for r in scored if r[field] == value}),
            )
            for value in sorted({r[field] for r in scored})
        }

    per_molecule = defaultdict(list)
    for row in scored:
        per_molecule[row["id"]].append(abs(row["error_V"]))
    macro_mae = (
        sum(sum(v) / len(v) for v in per_molecule.values()) / len(per_molecule)
        if per_molecule
        else None
    )
    spacing = []
    paired = defaultdict(list)
    for row in scored:
        paired[row["id"], row["doi"], row["source"]].append(row)
    for (gid, doi, source), rows in sorted(paired.items()):
        if len(rows) == 2 and {r["wave"] for r in rows} == {1, 2}:
            first, second = sorted(rows, key=lambda r: r["wave"])
            calc = first["calc_V"] - second["calc_V"]
            exp = first["exp_V"] - second["exp_V"]
            spacing.append(
                dict(
                    id=gid,
                    doi=doi,
                    source=source,
                    calc_V=calc,
                    exp_V=exp,
                    error_V=calc - exp,
                )
            )
    return dict(
        schema=1,
        panel_sha256=panel["panel_sha256"],
        status="complete_descriptive_benchmark"
        if not pending
        else "incomplete_benchmark",
        target_events=len(targets),
        admitted_events=len(scored),
        pending_events=len(pending),
        target_molecules=len({r["id"] for r in targets.values()}),
        admitted_molecules=len(per_molecule),
        original_studies=len({r["doi"] for r in scored}),
        metrics=metrics([r["error_V"] for r in scored]),
        molecule_macro_MAE_V=macro_mae,
        by_family=group("family"),
        by_wave=group("wave"),
        coverage_by_family={
            name: dict(
                target_events=sum(r["family"] == name for r in targets.values()),
                admitted_events=sum(r["family"] == name for r in scored),
            )
            for name in sorted({r["family"] for r in targets.values()})
        },
        source_counts=dict(Counter(r["doi"] for r in scored)),
        signed_wave_spacings=spacing,
        uncertainty=[
            cluster_intervals(
                scored,
                field,
                replicates=panel["bootstrap_replicates"],
                seed=panel["bootstrap_seed"],
                minimum=panel["minimum_clusters_for_interval"],
            )
            for field in ("id", "doi")
        ],
        outliers=sorted(scored, key=lambda r: abs(r["error_V"]), reverse=True),
        pending=pending,
        fitted_correction=False,
        held_out_from_model_development=False,
        exact_graft_accuracy_validated=False,
        limits=[
            "Parent measurements do not validate exact grafts or polymer interactions",
            (
                "Cluster intervals describe resampling of the observed panel, "
                "not prediction intervals"
            ),
            (
                "Molecule clusters may share an original study; "
                "study clustering is reported separately"
            ),
            "Missing calculations are excluded from metrics and retained in coverage",
            (
                "No empirical correction, outlier removal or acceptance threshold "
                "chosen from these errors"
            ),
        ],
    )


def review_report(panel, review):
    """Read a checked snapshot with its full inventory, rejecting concurrent changes."""
    validate_panel(panel)
    review = Path(review)
    inputs = json.loads((review / "input_inventory.json").read_text())
    verify_hashes(inputs)
    summary = json.loads((review / "summary.json").read_text())
    checks = json.loads((review / "checks.json").read_text())
    states = json.loads((review / "states.json").read_text())
    if (
        summary["status"] not in ("checked", "partial_checked")
        or summary["arithmetic_discrepancies"] != 0
        or summary["changed_inputs"]
        or not checks
        or any(c.get("ok") is not True for c in checks)
    ):
        raise ValueError("Independent benchmark review failed")
    if (
        set(summary["selected_ids"]) != {r["id"] for r in panel["targets"]}
        or summary["requested_tier_A_potentials"] != len(panel["targets"])
        or summary["arithmetic_checks"] != len(checks)
    ):
        raise ValueError("Independent review scope/counts disagree with frozen panel")
    paths = [
        review / name
        for name in (
            "summary.json",
            "checks.json",
            "comparisons.csv",
            "input_inventory.json",
            "states.json",
        )
    ]
    snapshot_hashes = {str(p.resolve()): digest(p) for p in paths}
    with (review / "comparisons.csv").open(newline="") as stream:
        report = summarize(panel, list(csv.DictReader(stream)))
    if report["admitted_events"] != summary["tier_A_potential_events"]:
        raise ValueError("Independent admission count changed")
    for row in report["outliers"]:
        ox, red = row["event"].split("->")
        left, right = (states[row["id"] + ":" + state] for state in (ox, red))
        if left["status"] != "ready" or right["status"] != "ready":
            raise ValueError("Endpoint not independently ready")
        expected = left["G_eV"] - right["G_eV"] - summary["native_Fc_absolute_V"]
        if abs(row["calc_V"] - expected) > 1e-9:
            raise ValueError(
                "Comparison differs from independently checked endpoint free energies"
            )
    verify_hashes(inputs)
    verify_hashes(snapshot_hashes)
    report["review_snapshot_sha256"] = snapshot_hashes
    return report
