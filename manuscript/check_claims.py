"""Independently reconstruct the numerical statements in the frozen paper snapshot."""

import csv
import hashlib
import json
import math
import statistics as descriptive
from collections import Counter
from pathlib import Path


def audit_quantum_deployment(load):
    """Reconstruct deployment comparisons without using a positive status alone."""
    try:
        report = load('quantum_environment.json')
    except FileNotFoundError:
        return None
    states = load('quantum_environment_states.json')['states']
    changes = load('quantum_environment_compatibility.json')['third_party_changes']
    if (report['status'] != 'fresh_compatible_GPU_deployment_independently_checked'
            or report['original_environment_modified'] or report['method_packages_changed']
            or report['full_portability_claimed'] or len(changes) != 14
            or report['compatible_CUDA_runtime_changes'] != len(changes)
            or report['libxtb_sha256'] != 'c405797f34ed96470f23158406e6d8a8e102681d3c3cbbaea16692f514f22853'
            or any(not r['name'].startswith('nvidia-') or not r['name'].endswith('-cu12')
                or r['recorded'] == r['replacement'] for r in changes)):
        raise ValueError('Deployment profile identity or qualification scope changed')
    expected = {'formaldehyde': 1, 'methyl_radical': 2}
    if len(states) != 2 or {r['case'] for r in states} != set(expected):
        raise ValueError('Both deployment spin controls are required')
    comparisons = {r['case']: r for r in report['controls']}
    for row in states:
        case = row['case']
        old, new = row['original'], row['fresh']
        if (old['mult'] != expected[case] or new['mult'] != expected[case]
                or old['symbols'] != new['symbols']
                or sorted(old['symbols']) != sorted(['C','H','H','O'] if case == 'formaldehyde' else ['C','H','H','H'])
                or any(r['status'] != 'normal_primary_electronic_and_native_admission_passed' for r in (old, new))):
            raise ValueError('Deployment control identity/admission changed')
        for state in (old, new):
            if (len(state['coordinates_A']) != len(state['symbols'])
                    or any(len(x) != 3 or not all(math.isfinite(v) for v in x) for x in state['coordinates_A'])
                    or len(state['thermal']['frequencies_cm']) != 3*len(state['symbols'])-6
                    or not all(math.isfinite(v) for v in state['thermal']['frequencies_cm'])):
                raise ValueError('Deployment geometry or mode inventory is incomplete')
            spin = -8.617333262145e-5 * 298.15 * math.log(state['mult'])
            G = state['energy']['e_smd_eV'] + state['thermal']['nuclear_by_cutoff_eV']['50'] + spin
            if not math.isfinite(G) or abs(G - state['G_eV']) > 1e-9:
                raise ValueError('Deployment free-energy reconstruction failed')
        errors = dict(energy_eV=max(abs(old['energy'][k]-new['energy'][k]) for k in ('e_smd_eV','e_gas_eV')),
            free_energy_eV=abs(old['G_eV']-new['G_eV']),
            frequency_cm=max(abs(a-b) for a,b in zip(sorted(old['thermal']['frequencies_cm']), sorted(new['thermal']['frequencies_cm']))))
        direct_rms = math.sqrt(sum(sum((a-b)**2 for a,b in zip(x,y))
            for x,y in zip(old['coordinates_A'],new['coordinates_A'])) / len(old['symbols']))
        if direct_rms > report['tolerances']['aligned_RMSD_Angstrom']:
            raise ValueError('Deployment direct-coordinate bound failed')
        if any(not math.isfinite(v) or v > report['tolerances'][k] for k,v in errors.items()):
            raise ValueError('Deployment comparison exceeds the declared bounds')
        if any(abs(v-comparisons[case]['differences'][k]) > 1e-9 for k,v in errors.items()):
            raise ValueError('Deployment reported comparison does not reconstruct')
    return dict(controls=2, external_distributions=report['external_distributions'],
        compatible_CUDA_runtime_changes=len(changes), universal_portability_claimed=False)


def audit_reviewer_package(load):
    """Check release coverage against its frozen inventory and replay receipts."""
    try:
        current = load("review_package.json")
    except FileNotFoundError:
        return None  # Earlier immutable epochs predate this packaging report.
    package = load("review_package_manifest.json")
    data = load("review_package_data.json")
    replay = load("review_package_replay.json")
    arithmetic = load("review_package_arithmetic.json")
    hessians = load("review_package_hessians.json")
    states = {r["id"]: r for r in data["states"]}
    expected_hessians = {k for k, v in states.items() if v["partition_kind"] != "monatomic"}
    actual_hessians = {r["id"] for r in hessians["states"]}
    electronic = data["auxiliary"]["electronic"]["records"]
    if (len(states) != len(data["states"]) or len(states) != package["raw_states"]
            or len(states) != current["raw_states"] or len(states) != arithmetic["raw_states"]
            or len(electronic) != package["electronic_only_records"]
            or "package_manifest.json" in package["files_sha256"]
            or len(package["files_sha256"]) + 1 != current["archive_files"]):
        raise ValueError("Reviewer package inventory counts differ")
    if (actual_hessians != expected_hessians or len(actual_hessians) != len(hessians["states"])
            or len(actual_hessians) != hessians["checked"]
            or len(actual_hessians) != current["CPU_Hessian_states"]):
        raise ValueError("Reviewer package Hessian coverage differs")
    if (replay["package_manifest_sha256"] != current["package_manifest_sha256"]
            or replay["arithmetic_interpreter"] != "/usr/bin/python3 -I -S"
            or not replay["original_checkout_and_calculation_opens_blocked_by_audit_hook"]
            or not replay["optional_Hessian_reconstruction"]
            or current["public_archive_deposited"] != package["public_archive_deposited"]
            or current["public_archive_deposited"] != replay["code_and_dataset_publicly_deposited"]
            or current["public_DOI_created"] or replay["full_GPU_environment_restored"]):
        raise ValueError("Reviewer package replay/public-deposition scope differs")
    if (arithmetic["original_cluster_paths_opened"] or arithmetic["quantum_calculations_run"]
            or arithmetic["independent_physical_QC_reestablished"] or arithmetic["experimental_accuracy_inferred"]
            or hessians["original_cluster_paths_opened"] or hessians["quantum_calculations_run"]
            or hessians["xTB_library_loaded"] or hessians["original_physical_QC_repeated"]
            or hessians["experimental_accuracy_inferred"]):
        raise ValueError("Cached replay was promoted to fresh physical or experimental validation")
    coverage = current["primary_coverage"]
    for key in ("benchmark_endpoints", "reference_reactions", "base_assisted_values",
                "relative_basicity_controls", "finite_graph_partitions"):
        if arithmetic[key] != coverage[key]:
            raise ValueError("Reviewer package primary replay coverage differs")
    for key in ("benchmark_endpoints", "reference_reactions", "base_assisted_values"):
        if package[key] != coverage[key]:
            raise ValueError("Reviewer package manuscript inventory differs")
    if arithmetic["auxiliary_scientific_records"] != current["auxiliary_coverage"]:
        raise ValueError("Reviewer package auxiliary replay coverage differs")
    negative = replay["negative_checks"]
    if len(negative) != current["negative_checks"] or len({r["case"] for r in negative}) != len(negative) or not all(r["rejected"] for r in negative):
        raise ValueError("Reviewer package negative checks differ")
    frequency_error = max(r["frequency_max_error_cm"] for r in hessians["states"])
    nuclear_error = max(abs(x) for r in hessians["states"] for x in r["nuclear_errors_eV"].values())
    if (frequency_error != current["maximum_frequency_reconstruction_error_cm"]
            or nuclear_error != current["maximum_nuclear_reconstruction_error_eV"]):
        raise ValueError("Reviewer package reconstruction errors differ")
    try:
        draft = load("review_package_draft.json")
    except FileNotFoundError:
        draft = None
    if draft is not None:
        assets = {r["name"]: r for r in draft["assets"]}
        archive = Path(current["compressed_archive"]).name
        if (draft["draft"] is not True or draft["published_at"] is not None
                or draft["public_raw_archive_deposited"] or draft["DOI_created"]
                or assets[archive]["sha256"] != current["archive_sha256"]
                or assets[archive]["bytes"] != current["compressed_bytes"]
                or assets["package_manifest.json"]["sha256"] != current["package_manifest_sha256"]):
            raise ValueError("Reviewer archive draft status or asset identity differs")
    return dict(raw_states=len(states), electronic_only_records=len(electronic),
                archive_files=len(package["files_sha256"]) + 1, Hessian_states=len(actual_hessians),
                negative_checks=len(negative), public_deposition=False, full_GPU_restoration=False)


def audit_thermal_settings(load):
    """Recompute published cutoff aggregates; retain the primary 50 cm^-1 values."""
    try:
        result = load("thermal_cutoff_sensitivity.json")
    except FileNotFoundError:
        return None
    if (result["quantum_calculations_run"] or result["primary_results_modified"]
            or result["experimental_accuracy_inferred"] or result["candidate_prediction_intervals_established"]
            or result["settings"]["cutoffs_cm"] != [25, 50, 100]
            or result["settings"]["primary_cutoff_cm"] != 50):
        raise ValueError("Cutoff scenarios were promoted beyond thermal-setting sensitivity")
    rows = result["benchmark"]["endpoints"]
    original = {(r["id"], r["event"]): r for r in load("benchmark.json")["outliers"]}
    if len(rows) != len(original) or {(r["id"], r["event"]) for r in rows} != set(original):
        raise ValueError("Cutoff benchmark coverage or identity differs")
    for row in rows:
        source = original[(row["id"], row["event"])]
        if (row["exp_V"] != source["exp_V"] or abs(row["E_by_cutoff_V"]["50"]-source["calc_V"]) > 1e-9
                or row["wave"] != (1 if "neu" in row["event"].split("->") else 2)):
            raise ValueError("Primary benchmark or wave labels changed in cutoff analysis")
        span = 1000*(max(row["E_by_cutoff_V"].values())-min(row["E_by_cutoff_V"].values()))
        if abs(span-row["range_mV"]) > 1e-9:
            raise ValueError("Endpoint cutoff span differs")
    for key in ("25", "50", "100"):
        errors = [r["E_by_cutoff_V"][key]-r["exp_V"] for r in rows]
        computed = dict(MAE_V=math.fsum(abs(v) for v in errors)/len(errors),
            RMSE_V=math.sqrt(math.fsum(v*v for v in errors)/len(errors)), bias_V=math.fsum(errors)/len(errors))
        if any(abs(value-result["benchmark"]["metrics_by_cutoff"][key][name]) > 1e-12 for name,value in computed.items()):
            raise ValueError("Cutoff benchmark statistics do not reconstruct")
    largest = max(r["range_mV"] for r in rows)
    if abs(largest-result["benchmark"]["max_endpoint_range_mV"]) > 1e-9:
        raise ValueError("Maximum cutoff benchmark span differs")
    reactions = result["grafting"]["reactions"]
    if len({r["name"] for r in reactions}) != 87 or len(reactions) != 87:
        raise ValueError("Expected the 67 reference and 20 base-assisted cutoff reactions")
    for row in reactions:
        values = row["SI_deltaG_by_cutoff_kJ_mol"].values()
        if abs(max(values)-min(values)-row["range_kJ_mol"]) > 1e-9:
            raise ValueError("Reaction cutoff span differs")
    crossing = [r["name"] for r in reactions if min(r["SI_deltaG_by_cutoff_kJ_mol"].values()) < 0 < max(r["SI_deltaG_by_cutoff_kJ_mol"].values())]
    max_reaction = max(r["range_kJ_mol"] for r in reactions)
    if crossing != result["grafting"]["zero_crossings"] or abs(max_reaction-result["grafting"]["max_range_kJ_mol"]) > 1e-9:
        raise ValueError("Cutoff reaction sign/range inventory differs")
    review = load("thermal_cutoff_review.json")
    if (review["status"] != "thermal_cutoff_full_state_reconstruction_checked"
            or review["potential_sums"] != 3*len(rows) or review["reaction_sums"] != 3*len(reactions)
            or review["maximum_potential_reconstruction_error_V"] > 1e-8
            or review["maximum_reaction_reconstruction_error_kJ_mol"] > 1e-6):
        raise ValueError("Independent thermal-cutoff source review differs")
    return dict(benchmark_endpoints=len(rows), reaction_models=len(reactions),
                maximum_endpoint_span_mV=largest, maximum_reaction_span_kJ_mol=max_reaction,
                reaction_sign_crossings=len(crossing), prediction_intervals_established=False)


def audit_decision_followup(load):
    """Recompute finite partitions, capacity bounds and decision counts independently."""
    baseline = load("finite_ensemble.json")
    states = {r["key"]: r for r in baseline["states"]}
    if len(states) != 12 or baseline["checked_native_G"] != 12 or baseline["pending_keys"]:
        raise ValueError("Incomplete twelve-state finite baseline")
    conversion = 1.602176634e-19 * 6.02214076e23 / 1000
    rt = 8.31446261815324e-3 * 298.15 / conversion
    partitions = {}
    for pair in baseline["finite_pairs"]:
        result = {}
        for cutoff in ("25", "50", "100"):
            energies = [states[k]["G_by_rotor_cutoff_eV"][cutoff] for k in pair["conformer_keys"]]
            minimum = min(energies)
            g = minimum - rt * math.log(math.fsum(math.exp(-(e-minimum)/rt) for e in energies))
            if abs(g-pair["independent_G_eV_by_rotor_cutoff_cm"][cutoff]) > 1e-8:
                raise ValueError("Finite conformer partition does not reconstruct")
            result[cutoff] = g
        partitions[pair["graph_key"]] = result
    for row in baseline["cycles"]:
        i, j = row["sites"]
        for cutoff in ("25", "50", "100"):
            dg = math.fsum((partitions[f"occ_{i}_{j}"][cutoff], partitions["occ_none"][cutoff],
                           -partitions[f"occ_{i}"][cutoff], -partitions[f"occ_{j}"][cutoff])) * conversion
            if abs(dg-row["deltaDeltaG_kJ_mol_by_rotor_cutoff_cm"][cutoff]) > 1e-6:
                raise ValueError("Finite ensemble occupancy cycle does not reconstruct")
        for terms in row["conditional_DBU_grafting"].values():
            for second, first in (("graft_j_given_i_kJ_mol", "graft_j_on_empty_kJ_mol"),
                                  ("graft_i_given_j_kJ_mol", "graft_i_on_empty_kJ_mol")):
                if abs(terms[second]-terms[first]-row["deltaDeltaG_kJ_mol_by_rotor_cutoff_cm"]["50"]) > 1e-6:
                    raise ValueError("Finite conditional grafting does not close")
    capacity = load("pmdi_capacity.json")
    if capacity["accessible_capacity_established"] or capacity["experimental_loading_inferred"]:
        raise ValueError("Composition bounds were promoted to material performance")
    weights = dict(C=12.011, H=1.008, N=14.007, O=15.999, Cl=35.453)
    faraday = conversion * 1000
    graph_bounds = {}
    for graph in capacity["exact_graphs_at_handle_density_0_6"]:
        grafts = len(graph["occupied_sites"])
        atoms = dict(C=43+12*grafts, H=45+7*grafts, Cl=3-grafts, N=2*grafts, O=4*grafts)
        formula = "".join(k+(str(v) if v != 1 else "") for k, v in atoms.items() if v)
        mass = math.fsum(weights[k]*v for k, v in atoms.items())
        if graph["formula"] != formula or abs(mass-graph["mass_g_mol"]) > 1e-8:
            raise ValueError("Exact occupancy mass does not reconstruct")
        for n, label in ((1, "one_e"), (2, "two_e")):
            expected = grafts*n*faraday/(3.6*mass)
            if abs(expected-graph[label+"_upper_bound_mAh_g"]) > 1e-8:
                raise ValueError("Capacity upper bound does not reconstruct")
        graph_bounds[tuple(graph["occupied_sites"])] = graph["two_e_upper_bound_mAh_g"]
    combinations = set()
    for row in capacity["composition_grid"]:
        n, h, g = (row[k] for k in ("repeat_units", "reactive_handles", "grafts"))
        if not 0 <= g <= h <= n or (n, h, g) in combinations:
            raise ValueError("Invalid composition grid")
        combinations.add((n, h, g))
        # Independently sum atomic compositions; do not import the producer's mass model.
        atoms = dict(C=8*n+h+12*g, H=8*n+h+2+7*g, Cl=h-g, N=2*g, O=4*g)
        mass = math.fsum(weights[k]*v for k, v in atoms.items())
        for electrons in (1, 2):
            result = row[f"{electrons}e"]
            expected = g*electrons*faraday/(3.6*mass)
            if (abs(result["mass_g_mol"]-mass) > 1e-8
                    or abs(result["capacity_upper_bound_mAh_g"]-expected) > 1e-8):
                raise ValueError("Composition capacity grid does not reconstruct")
    summary = load("all_axes_summary.json")
    decisions = load("all_axes_decisions.json")["candidates"]
    required = [r["id"] for r in load("all_axes_policy.json")["required_axes"]]
    if len({r["product_key"] for r in decisions}) != len(decisions) or len(decisions) != summary["exact_products"]:
        raise ValueError("Decision identities or counts differ")
    if dict(Counter(r["status"] for r in decisions)) != summary["decisions"]:
        raise ValueError("Decision status counts differ")
    eligible = []
    for row in decisions:
        axes = {r["axis"]: r["status"] for r in row["axes"]}
        if set(axes) != set(required) or len(row["axes"]) != len(required):
            raise ValueError("Required decision axes differ")
        passed = all(axes[k] == "pass" for k in required)
        if row["all_required_axes_pass"] != passed or row["eligible_for_Pareto"] != passed:
            raise ValueError("Candidate bypassed all-axis admission")
        if passed:
            eligible.append(row["id"])
    if sorted(eligible) != sorted(summary["eligible_ids"]):
        raise ValueError("Eligible candidate counts differ")
    precursor = load("pmdi_precursor_review.json")
    route = precursor["prospective_route_decision"]
    if (not route["evidence"]["exact_precursor_access_reviewed"] or route["route_status"] != "hold"
            or route["proceed_to_full_redox"] or precursor["exact_polymer_grafting_demonstrated"]):
        raise ValueError("Precursor access was confused with a completed attachment route")
    deformation = load("segment_relaxation.json")
    for pair in deformation["pairs"]:
        if abs(math.fsum(f["deformation_kJ_mol"] for f in pair["fragments"])
               - pair["deformation_kJ_mol_CODATA"]) > 1e-6:
            raise ValueError("Signed fragment deformation does not close")
    origin = None
    try:
        origin = load("basicity_reference_origin.json")
    except FileNotFoundError:
        pass  # Earlier immutable epochs predate this source-origin review.
    if origin is not None:
        if (origin["independent_experimental_pKaH_MeCN"] is not None
                or origin["source_role"]["independent_experimental_origin_eligible"]
                or origin["source"]["origin_class"] != "computed_PCM_estimate"
                or origin["primary_grafting_results_modified"] or origin["candidate_decisions_modified"]):
            raise ValueError("Computed DIPEA source was promoted to experimental validation")
        capture = -.00831446261815324*298.15*math.log(10)*(origin["source"]["reported_value"]-load("acid_base.json")["acid_anchor"]["pKa_MeCN"])
        if abs(capture-origin["computed_literature_capture_kJ_mol"]) > 1e-10:
            raise ValueError("DIPEA model sensitivity capture does not reconstruct")
        for row in origin["candidate_model_sensitivities"]:
            if abs(row["reference_deltaG_kJ_mol"]+capture-row["fixed_reference_plus_computed_literature_capture_kJ_mol"]) > 1e-9 or row["calibrated_full_route_prediction"]:
                raise ValueError("DIPEA sensitivity changed or acquired a calibration claim")
    return dict(finite_native_G=len(states), finite_graph_partitions=len(partitions),
                finite_cycles=len(baseline["cycles"]), capacity_graphs=len(graph_bounds),
                capacity_compositions=len(combinations),
                all_axis_graphs=len(decisions), all_axis_eligible=len(eligible),
                all_axis_status_counts=summary["decisions"], signed_deformation_pairs=len(deformation["pairs"]))


def audit():
    root = Path(__file__).resolve().parents[1]
    pointer = json.loads((root / "manuscript/evidence/current.json").read_text())
    folder = (root / pointer["manifest"]).parent
    load = lambda name: json.loads((folder / name).read_text())
    benchmark = load("benchmark.json")
    panel = load("benchmark_panel.json")
    rows = list(csv.DictReader((folder / "benchmark_points.csv").open()))
    errors = [float(r["calc_V"]) - float(r["exp_V"]) for r in rows]
    keys = [(r["id"], r["event"], r["doi"]) for r in rows]
    targets = {(r["id"], r["event"], r["doi"]): r for r in panel["targets"]}
    if len(set(keys)) != len(keys) or not set(keys) <= set(targets):
        raise ValueError("Duplicated or non-panel manuscript benchmark point")
    if any(abs(float(r["exp_V"]) - targets[key]["exp_V"]) > 1e-12 for r, key in zip(rows, keys)):
        raise ValueError("Manuscript experimental label differs from frozen panel")
    statistics = dict(MAE_V=math.fsum(abs(e) for e in errors) / len(errors),
                      RMSE_V=math.sqrt(math.fsum(e * e for e in errors) / len(errors)),
                      bias_V=math.fsum(errors) / len(errors), maximum_absolute_error_V=max(abs(e) for e in errors))
    if any(abs(statistics[key] - benchmark["metrics"][key]) > 1e-12 for key in statistics):
        raise ValueError("Manuscript benchmark arithmetic mismatch")
    if (len(rows) != benchmark["admitted_events"] or len(targets) != benchmark["target_events"]
            or len({r["id"] for r in rows}) != benchmark["admitted_molecules"]
            or len({r["doi"] for r in rows}) != benchmark["original_studies"]):
        raise ValueError("Manuscript coverage mismatch")
    bias = load("reference_bias.json")
    if bias["primary_potentials_modified"] or bias["external_validation_claimed"]:
        raise ValueError("Diagnostic offsets replaced the primary benchmark")
    cv_checks = 0
    for model in bias["cross_validation"]:
        residuals = []
        for index, row in enumerate(rows):
            label = model["holdout"]
            members = [j for j, other in enumerate(rows)
                       if (j != index if label == "endpoint" else other[label] != row[label])
                       and (model["model"] == "constant_offset"
                            or targets[keys[j]]["wave"] == targets[keys[index]]["wave"])]
            if not members:
                continue
            offset = math.fsum(errors[j] for j in members) / len(members)
            residuals.append(errors[index] - offset)
        checked = dict(MAE_V=math.fsum(abs(e) for e in residuals)/len(residuals),
                       RMSE_V=math.sqrt(math.fsum(e*e for e in residuals)/len(residuals)),
                       bias_V=math.fsum(residuals)/len(residuals))
        if any(abs(checked[k]-model["metrics"][k]) > 1e-12 for k in checked):
            raise ValueError("Grouped offset diagnostic does not reconstruct")
        cv_checks += 1
    for interval in bias["uncertainty"]:
        primary = next(r for r in benchmark["uncertainty"] if r["cluster"] == interval["cluster"])
        if interval["intervals"] != primary["intervals"]:
            raise ValueError("Printed bootstrap interval differs from primary evidence")
    base = load("base_assisted.json")
    base_values = 0
    for candidate in base["candidates"]:
        for name, modes in candidate["base_models"].items():
            for mode, value in modes.items():
                if value is None:
                    continue
                if abs(value - candidate["reference_deltaG_kJ_mol"] - base["corrections"][name][mode]) > 1e-7:
                    raise ValueError("Manuscript base-assisted cycle mismatch")
                base_values += 1
    polymer = load("polymer_validation.json")
    cases = {r["key"]: r for r in polymer["cases"]}
    conversion_differences = []
    for cycle in polymer["cycles"]:
        difference_eV = math.fsum(cases[key]["E_eV"] * weight for key, weight in cycle["coefficients"].items())
        # Existing source-frozen records use this documented five-digit factor.
        # Also record its difference from the SI conversion rather than hiding it.
        value = difference_eV * 96.485
        conversion_differences.append(abs(difference_eV * 96.48533212 - value))
        if abs(value - cycle["electronic_deltaDeltaE_kJ_mol"]) > 1e-6:
            raise ValueError("Manuscript electronic occupancy cycle mismatch")
        if any(cases[key]["G_eV"] is None for key in cycle["coefficients"]):
            if cycle["free_energy_deltaDeltaG_kJ_mol"] is not None:
                raise ValueError("Missing native thermochemistry was filled in")
        else:
            dg = math.fsum(cases[key]["G_eV"] * weight for key, weight in cycle["coefficients"].items()) * 96.485
            if abs(dg - cycle["free_energy_deltaDeltaG_kJ_mol"]) > 1e-6:
                raise ValueError("Manuscript native occupancy cycle mismatch")
            for cutoff, expected in cycle["rotor_cutoff_sensitivity_kJ_mol"].items():
                value = math.fsum(cases[key]["cutoff_G_eV"][cutoff] * weight for key, weight in cycle["coefficients"].items()) * 96.485
                if abs(value - expected) > 1e-6:
                    raise ValueError("Manuscript rotor-cutoff sensitivity mismatch")
            for terms in cycle["conditional_grafting_kJ_mol"].values():
                for second, first in (("second_j_given_i", "first_j"), ("second_i_given_j", "first_i")):
                    if abs(terms[second] - terms[first] - dg) > 1e-6:
                        raise ValueError("Conditional grafting difference does not close")
    routes = load("exact_routes.json")
    if routes["completed"] != routes["expected"] or sum(routes["counts"].values()) != routes["completed"]:
        raise ValueError("Exact route handoff coverage mismatch")
    acid = load("acid_base.json")
    for row in acid["rows"]:
        if row["reference"] is None:
            if row["external_capture_kJ_mol"] is not None:
                raise ValueError("Missing exact-base anchor filled in")
            continue
        expected = -.00831446261815324 * acid["temperature_K"] * math.log(10) * (row["reference"]["pKaH_MeCN"] - acid["acid_anchor"]["pKa_MeCN"])
        if abs(expected - row["external_capture_kJ_mol"]) > 1e-10:
            raise ValueError("External acid-capture conversion mismatch")
    concurrency = load("cpu_concurrency.json")
    if concurrency["run_count"] != len(concurrency["runs"]) or concurrency["serial_discrepancies"]:
        raise ValueError("Concurrency experiment count or independent comparison differs")
    for condition in concurrency["conditions"]:
        runs = [r for r in concurrency["runs"] if r["condition"] == condition["condition"]]
        if len(runs) != condition["repeats"]:
            raise ValueError("Concurrency repeat count differs")
        for key in ("launch_to_observed_completion_seconds", "workers_executing_chunks"):
            values = [r[key] for r in runs]
            expected = dict(median=descriptive.median(values), minimum=min(values), maximum=max(values))
            if expected != condition[key]:
                raise ValueError("Concurrency summary does not reconstruct from actual runs")
    fresh = load("fresh_checkout.json")
    if not fresh["figure_CSV_exact_match"] or not fresh["figure_RGB_pixels_exact_match"] or fresh["serial_discrepancies"]:
        raise ValueError("Fresh-checkout example did not reproduce the checked results")
    transfer = load("proton_transfer.json")
    components = {r["id"]: r["G_eV"] for r in transfer["states"]}
    references = load("proton_transfer_sources.json")
    for row, anchor in zip(transfer["comparisons"], references["comparisons"]):
        if row["base"] != anchor["base"]:
            raise ValueError("Relative-basicity source identity mismatch")
        key = row["base"]
        dg = (components[key+"_protonated"] + components["dbu_neutral"]
              - components[key+"_neutral"] - components["dbu_protonated"]) * 96.48533212331002
        external = .00831446261815324 * 298.15 * math.log(10) * (
            anchor["reference"]["pKaH_MeCN"]-anchor["source"]["pKaH_MeCN"])
        if abs(dg-row["DFT_deltaG_kJ_mol"]) > 1e-6 or abs(external-row["external_deltaG_kJ_mol"]) > 1e-9:
            raise ValueError("Relative-basicity comparison does not reconstruct")
    restart = load("quantum_restart.json")
    for control in restart["controls"]:
        if any(control["differences"][key] > tolerance for key, tolerance in restart["tolerances"].items()):
            raise ValueError("Quantum restart exceeds its predeclared tolerance")
    basis = load("neutral_basis.json")
    if basis["checked"] != 5 or basis["pending"] or basis["primary_values_modified"]:
        raise ValueError("Neutral basis panel is incomplete or modified primary values")
    for row in basis["rows"]:
        if abs(row["delta_thermal_eV"] + row["delta_geometry_SP_eV"]
               - row["delta_first_reduction_V_holding_anion_and_Fc_fixed"]) > 1e-8:
            raise ValueError("Basis sensitivity components do not close")
    exact = load("exact_case_join.json")
    for row in exact["cases"]:
        g = row["state_G_eV"]
        if (abs(g["neu"] - g["red1"] - exact["Fc_absolute_V"] - row["E1_V_vs_Fc"]) > 1e-9
                or abs(g["red1"] - g["red2"] - exact["Fc_absolute_V"] - row["E2_V_vs_Fc"]) > 1e-9
                or row["route"]["proceed_to_full_redox"] or row["end_to_end_complete"]):
            raise ValueError("Exact-case values or pending route status changed")
    cached = load("cached_states_manifest.json")
    states = {row["id"]: row["checked_G_eV"] for row in cached["states"]}
    reproduced = load("cached_reproduction.json")
    p = cached["potential"]
    potential = states[p["oxidized"]] - states[p["reduced"]] - states[p["reference_oxidized"]] + states[p["reference_reduced"]]
    if abs(potential - reproduced["PMDI_first_reduction_V_vs_Fc"]) > 1e-9:
        raise ValueError("Cached example potential changed")
    for recipe, value in zip(cached["reactions"], reproduced["reactions"]):
        dg = math.fsum(states[key] * nu for key, nu in recipe["coefficients"].items())
        if abs(dg - value["deltaG_eV"]) > 1e-9:
            raise ValueError("Cached example grafting sum changed")
    added = audit_decision_followup(load) if (folder / "finite_ensemble.json").exists() else None
    release = audit_reviewer_package(load)
    cutoff = audit_thermal_settings(load)
    deployment = audit_quantum_deployment(load)
    if cutoff is not None:
        if hashlib.sha256((folder / "thermal_cutoff_sensitivity.json").read_bytes()).hexdigest() != load("thermal_cutoff_review.json")["result_sha256"]:
            raise ValueError("Frozen thermal-cutoff result differs from its independent source review")
    if release is not None:
        package = load("review_package_manifest.json")
        current = load("review_package.json")
        digest = lambda name: hashlib.sha256((folder / name).read_bytes()).hexdigest()
        if (digest("review_package_manifest.json") != current["package_manifest_sha256"]
                or digest("review_package_data.json") != package["files_sha256"]["data/manifest.json"]
                or package["code_commit"] != current["code_snapshot_commit"]):
            raise ValueError("Frozen reviewer inventory differs from the tested package")
    return dict(status="numerical_snapshot_reconstructed", benchmark=dict(admitted=len(rows), target=len(targets),
                molecules=len({r["id"] for r in rows}), studies=len({r["doi"] for r in rows}), **statistics),
                base_assisted_values=base_values, qualified_spacing_G=sum(r["G_eV"] is not None for r in cases.values()),
                spacing_graphs=len(cases), exact_route_entries=len(load("route_evidence.json")["candidates"]),
                exact_routes_checked=routes["completed"], route_decisions=routes["counts"],
                acid_base_scope=acid["status"], conditional_native_cycles=len(polymer["cycles"]),
                checked_concurrency_runs=len(concurrency["runs"]), fresh_checkout_designs=fresh["exact_designs"],
                source_energy_conversion_kJ_mol_per_eV=96.485,
                maximum_SI_conversion_difference_kJ_mol=max(conversion_differences),
                checked_offset_models=cv_checks, relative_basicity_controls=len(transfer["comparisons"]),
                quantum_restart_controls=len(restart["controls"]),
                neutral_basis_controls=basis["checked"], exact_cases_with_explicit_holds=len(exact["cases"]),
                decision_followup=added,
                reviewer_package=release,
                thermal_cutoff_sensitivity=cutoff,
                quantum_environment=deployment,
                experimental_accuracy_inferred=False, submission_ready=False)


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2))
