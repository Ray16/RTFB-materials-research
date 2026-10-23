"""Conjunctive, evidence-bound candidate admission before Pareto comparison.

An acceptable value on another axis never rescues an unacceptable or unresolved
required axis. Bounds describe the declared model/scenarios, not a universal
experimental confidence interval. This module does not fill missing evidence or
choose application cutoffs from the distribution of candidates being screened.
"""
from __future__ import annotations

import hashlib
import json
import math
import re


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def policy_hash(policy):
    return hashlib.sha256(json.dumps(policy, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False).encode()).hexdigest()


def validate_policy(policy):
    axes = policy["required_axes"]
    if not policy.get("profile_id") or not axes:
        raise ValueError("A named profile and at least one required axis are needed")
    names = [axis["id"] for axis in axes]
    if len(set(names)) != len(names):
        raise ValueError("Required axes must be unique")
    for axis in axes:
        if not axis.get("scope") or not axis.get("criterion_source"):
            raise ValueError("Each criterion needs an explicit scope and source")
        kind = axis["kind"]
        if kind == "numeric":
            if not axis.get("unit"):
                raise ValueError("Numeric criteria need units")
            lo, hi = axis.get("minimum"), axis.get("maximum")
            if any(v is not None and not finite(v) for v in (lo, hi)):
                raise ValueError("Cutoffs must be finite numbers or explicitly unset")
            if lo is not None and hi is not None and lo > hi:
                raise ValueError("Acceptance range is reversed")
        elif kind == "categorical":
            accepted = axis.get("accepted_values", [])
            rejected = axis.get("rejected_values", [])
            if not accepted or set(accepted) & set(rejected):
                raise ValueError("Categorical acceptance and rejection must be distinct")
        else:
            raise ValueError("Unknown criterion kind")
    for group in policy.get("shared_context_groups", []):
        if not group.get("id") or not set(group.get("axes", [])) <= set(names) or not group.get("axes"):
            raise ValueError("Shared-context groups must name required axes")


def assess_axis(axis, record, product_key):
    """Pass only qualified, exact-identity evidence wholly within its criterion."""
    result = dict(axis=axis["id"], status="pending", reason="Required evidence missing")
    if not record:
        return result
    if record.get("product_key") != product_key:
        raise ValueError("Evidence belongs to a different exact product")
    if record.get("scope") != axis["scope"]:
        raise ValueError("Evidence has a different molecular/model scope")
    result["evidence_status"] = record.get("status", "missing")
    if record.get("status") != "qualified":
        result["reason"] = record.get("reason") or "Evidence is incomplete or unvalidated"
        return result
    sources = record.get("sources")
    if not isinstance(sources, list) or not sources or any(
        not isinstance(source, dict) or not source.get("path")
        or not isinstance(source.get("sha256"), str)
        or re.fullmatch(r"[0-9a-f]{64}", source["sha256"]) is None
        for source in sources
    ):
        raise ValueError("Qualified evidence needs source hashes")
    if axis["kind"] == "categorical":
        value = record.get("value")
        if value in axis["accepted_values"]:
            result.update(status="pass", reason="Reviewed evidence meets the criterion")
        elif value in axis.get("rejected_values", []):
            result.update(status="fail", reason="Reviewed evidence violates the criterion")
        else:
            result["reason"] = "Reviewed outcome remains unresolved for this criterion"
        return result
    if record.get("unit") != axis["unit"]:
        raise ValueError("Evidence and criterion units differ")
    lo, hi = axis.get("minimum"), axis.get("maximum")
    if lo is None and hi is None:
        result["reason"] = "Application cutoff has not been justified and declared"
        return result
    interval = record.get("interval")
    if record.get("interval_kind") not in axis.get("accepted_interval_kinds", []):
        result["reason"] = "No accepted uncertainty/sensitivity treatment for this axis"
        return result
    if not isinstance(interval, (list, tuple)) or len(interval) != 2:
        result["reason"] = "Complete uncertainty/sensitivity interval missing"
        return result
    a, b = interval
    if not finite(a) or not finite(b) or a > b:
        raise ValueError("Invalid evidence interval")
    value = record.get("value")
    if value is not None and (not finite(value) or not a <= value <= b):
        raise ValueError("Point estimate does not belong to its declared interval")
    result["interval"] = list(interval)
    if (lo is None or a >= lo) and (hi is None or b <= hi):
        result.update(status="pass", reason="Entire declared interval meets the criterion")
    elif (lo is not None and b < lo) or (hi is not None and a > hi):
        result.update(status="fail", reason="Entire declared interval violates the criterion")
    else:
        result["reason"] = "Interval crosses an acceptance boundary"
    return result


def assess_candidate(candidate, policy):
    validate_policy(policy)
    key = candidate.get("product_key")
    if not key:
        raise ValueError("An exact product identity is required")
    if candidate.get("role", "candidate") not in ("candidate", "comparison_control"):
        raise ValueError("Unknown candidate/control role")
    records = candidate.get("axes", {})
    unknown = set(records) - {a["id"] for a in policy["required_axes"]}
    if unknown:
        raise ValueError("Evidence contains axes outside the declared profile")
    axes = [assess_axis(axis, records.get(axis["id"]), key)
            for axis in policy["required_axes"]]
    by_axis = {a["axis"]: a for a in axes}
    for group in policy.get("shared_context_groups", []):
        contexts = []
        for name in group["axes"]:
            record = records.get(name) or {}
            if record.get("status") != "qualified":
                continue
            context = record.get("context")
            if not isinstance(context, dict) or not context:
                if by_axis[name]["status"] in ("pass", "fail"):
                    by_axis[name].update(status="pending", reason="Required operating/reaction context is not declared")
                continue
            for field, expected in group.get("required_context", {}).items():
                if field not in context:
                    by_axis[name].update(status="pending", reason="Required context field missing: " + field)
                elif context[field] != expected:
                    raise ValueError("Evidence belongs to a different operating/reaction context")
            contexts.append(context)
        if contexts and any(context != contexts[0] for context in contexts[1:]):
            raise ValueError("Required axes were evaluated in incompatible contexts")
    failed = [row["axis"] for row in axes if row["status"] == "fail"]
    pending = [row["axis"] for row in axes if row["status"] == "pending"]
    eligible = not failed and not pending and candidate.get("role", "candidate") == "candidate"
    status = "eligible" if eligible else "excluded" if failed else "pending"
    if candidate.get("role", "candidate") != "candidate":
        status = "comparison_control"
    return dict(id=candidate["id"], product_key=key, profile_id=policy["profile_id"],
                policy_sha256=policy_hash(policy), required_axis_ids=[a["id"] for a in policy["required_axes"]],
                status=status, all_required_axes_pass=not failed and not pending,
                eligible_for_Pareto=eligible, failed_axes=failed, pending_axes=pending,
                axes=axes, inference_scope=policy["decision_scope"],
                experimental_performance_established=False)


def feasible_ids(decisions, *, policy):
    """A front consumer must obtain eligibility here, not from descriptor rank."""
    validate_policy(policy)
    expected_hash = policy_hash(policy)
    expected_axes = {a["id"] for a in policy["required_axes"]}
    ids = set()
    for result in decisions:
        if (result.get("policy_sha256") != expected_hash
                or result.get("profile_id") != policy["profile_id"]
                or set(result.get("required_axis_ids", [])) != expected_axes
                or {a["axis"] for a in result.get("axes", [])} != expected_axes
                or len(result.get("axes", [])) != len(expected_axes)):
            raise ValueError("Decision belongs to another or incomplete criterion profile")
        if result.get("eligible_for_Pareto"):
            if (result.get("status") != "eligible"
                    or result.get("all_required_axes_pass") is not True
                    or not result.get("axes")
                    or any(axis.get("status") != "pass" for axis in result["axes"])):
                raise ValueError("Inconsistent full-axis admission record")
            if result["id"] in ids:
                raise ValueError("Duplicate decision identity")
            ids.add(result["id"])
    return ids


def pareto_after_feasibility(candidates, policy, objectives, *, pool=None):
    """Conservative interval dominance using the very same qualified evidence.

    A is no worse only when its least favorable bound is at least as good as B's
    most favorable bound on every objective. Overlapping uncertainty therefore
    cannot manufacture a confident preference. No scalarized score is used.
    """
    validate_policy(policy)
    if policy.get("context", {}).get("pools_evaluated_separately"):
        if pool not in ("anolyte", "catholyte"):
            raise ValueError("Select the electrode pool before comparing candidates")
        if any("electrode_pools" not in c for c in candidates):
            raise ValueError("Candidate electrode roles are not declared")
        candidates = [c for c in candidates if pool in c["electrode_pools"]]
    if len(objectives) < 2 or len({o["axis"] for o in objectives}) != len(objectives):
        raise ValueError("Declare at least two distinct Pareto objectives")
    rules = {a["id"]: a for a in policy["required_axes"]}
    for objective in objectives:
        if (objective["axis"] not in rules or rules[objective["axis"]]["kind"] != "numeric"
                or objective["direction"] not in ("minimize", "maximize", "minimize_magnitude")):
            raise ValueError("Pareto objectives must be declared numeric axes with a direction")
    if len({c["id"] for c in candidates}) != len(candidates):
        raise ValueError("Duplicate candidate identity in a comparison pool")
    decisions = [assess_candidate(candidate, policy) for candidate in candidates]
    eligible = feasible_ids(decisions, policy=policy)
    values = {}
    for candidate in candidates:
        if candidate["id"] not in eligible:
            continue
        values[candidate["id"]] = []
        for objective in objectives:
            record = candidate["axes"][objective["axis"]]
            a, b = record["interval"]
            if objective["direction"] == "minimize_magnitude":
                a, b = (0 if a <= 0 <= b else min(abs(a), abs(b))), max(abs(a), abs(b))
            values[candidate["id"]].append((a, b) if objective["direction"] == "maximize" else (-b, -a))
    def dominates(left, right):
        return (all(a[0] >= b[1] for a, b in zip(left, right))
                and any(a[0] > b[1] for a, b in zip(left, right)))
    dominated_by = {key: sorted(other for other in values if other != key
                               and dominates(values[other], values[key])) for key in values}
    return dict(profile_id=policy["profile_id"], policy_sha256=policy_hash(policy), electrode_pool=pool,
                eligible_ids=sorted(eligible), front_ids=sorted(k for k, v in dominated_by.items() if not v),
                dominated_by=dominated_by, objectives=objectives, decisions=decisions,
                inference_scope=policy["decision_scope"],
                rule="All-axis admission, then conservative interval dominance; no weighted score")
