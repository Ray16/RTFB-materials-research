"""Refresh navigation metadata from authoritative artifacts; preserve the old index."""
import hashlib,json,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
base=Path(__file__).resolve().parent
index=ROOT/'results/LATEST_REVIEW.json'
archive=base/'LATEST_REVIEW_before_sync.json'
if not archive.exists():archive.write_bytes(index.read_bytes())
feed_path=ROOT/'results/candidates/current.json';feed=json.loads(feed_path.read_text());s=feed['summary']
benchmark_path=ROOT/'results/reviews/2026-10-18/matched_parent_projected_checked/summary.json';benchmark=json.loads(benchmark_path.read_text())
global_path=base/'atomic_repair_adoptions_tables/review_status.json';global_review=json.loads(global_path.read_text());publication_path=global_path.parent/'publication.json';publication=json.loads(publication_path.read_text()) if publication_path.exists() else {}
slides_path=ROOT/'results/reviews/2026-10-18/slides_projected_publisher/status.json';slides=json.loads(slides_path.read_text())
thermal_path=ROOT/'results/validation/rigid_mode_thermal_2026-10-18/kernel_isolated_independent_review/summary.json';thermal=json.loads(thermal_path.read_text())
geometry_path=ROOT/'results/validation/geometry_driver_full_response_2026-10-18/status.json';geometry=json.loads(geometry_path.read_text())
assert s['arithmetic_discrepancies']==0 and thermal['failed_checks']==0 and benchmark['arithmetic_discrepancies']==0
current=dict(as_of_Chicago=feed['as_of_Chicago'],scientific_validation_complete=False,
 candidate_results=dict(guide='results/candidates/README.md',workbook='results/candidates/latest.xlsx',additional_csv='results/candidates/latest.csv',matched_comparisons='results/candidates/comparisons.csv',snapshot=feed['snapshot_path'],additional_candidates=s['additional_candidate_count'],original39_with_native_potentials=s['original39_native_potential_molecules'],complete_potential_pairs=s['paired_complete_potential_count'],arithmetic_checks=s['arithmetic_checks'],arithmetic_discrepancies=0,publication_status='corrected internal-mode thermochemistry; field-specific QC',inputs_current_during_collection=not s['inputs_changed_during_snapshot']),
 thermal_projection_review=dict(guide='results/validation/rigid_mode_thermal_2026-10-18/README.md',status='corrected_protocol_integrated_and_independently_admitted',qualified_cases=thermal['states'],independent_checks=thermal['checks'],independent_failed_checks=0,producer_and_both_primary_readers_integrated=True,soft_internal_modes_retained=True),
 native_Fc_absolute_V=benchmark['native_Fc_absolute_V'],
 matched_parent_validation=dict(path=str(benchmark_path.relative_to(ROOT)),scope=benchmark['scope'],potentials=benchmark['tier_A_potential_events'],molecules=benchmark['tier_A_potential_molecules'],studies=benchmark['tier_A_original_studies'],stats=benchmark['tier_A_potential_stats'],independent_checks=benchmark['arithmetic_checks'],exact_graft_accuracy_validated=False),
 geometry_driver=dict(path=str(geometry_path.relative_to(ROOT)),status=geometry['status'],full_grid_response=True,explicit_optimizer_and_SCF_convergence=True,actual_CPU_smoke=geometry['actual_cpu_smoke']),
 source_frozen_refinement=dict(plan='results/validation/ndi_methyl_radical_minimum_2026-10-19/source_plus.json',launch='results/validation/ndi_methyl_radical_minimum_2026-10-19/launch.json',actual_start_verification='results/validation/ndi_methyl_radical_minimum_2026-10-19/start_verification.json',rule='Do not edit the plan-bound numerical producers until both future derivative qualification checks have passed or their workers are verified terminal.'),
 radical_adoptions=dict(ndi_methoxyethyl='results/validation/ndi_radical_adoptions_2026-10-19/ndi_methoxyethyl/adoption_receipt.json',ndi_butyl='results/validation/ndi_radical_adoptions_2026-10-19/ndi_butyl/adoption_receipt.json',independent_checks_each=157),
 global_review=dict(path=str(global_path.relative_to(ROOT)),status=global_review['status'],arithmetic_checks=publication.get('arithmetic_checks'),arithmetic_discrepancies=global_review['arithmetic_discrepancies'],inputs_current_during_review=global_review.get('inputs_current'),changed_inputs=global_review.get('changed_inputs',[]),scientific_accuracy_validated=False),
 slides=dict(status=slides['status'],path='slides/redox_workflow.pdf',publisher=str(slides_path.relative_to(ROOT)),snapshot=slides.get('snapshot'),pdf_sha256=slides.get('pdf_sha256'),pages=slides.get('pages')),
 software_tests_passed=395,software_checks='results/reviews/2026-10-19/software_checks_atomic_repair.log',
 remaining_work=['Hexyl NDI native minimum/thermochemistry','CN-PMDI near-spin ensemble admission','Methyl NDI radical minimum refinement','Solution reorganization and kinetic qualification','Broader experimental benchmark and exact-graft validation'])
report=dict(schema=2,current=current,scientific_validation_complete=False,updated=time.time(),historical_metadata=str(archive.relative_to(ROOT)),progress_document='docs/SCIENTIFIC_PROGRESS_2026-10-18.md',discovery_source_catalog='config/discovery_sources.json',candidate_search='results/discovery/systematic_families_2026-10-18/summary.json',candidate_results_entry='results/candidates/README.md',thermal_publication_review='config/thermal_publication_review.json',note='Current scoped evidence is separate from broad scientific validation. Historical unprojected references and counts are preserved in historical_metadata; no parent measurement is assigned to an exact graft.')
atomic_write_json(index,report)
quality=dict(status='corrected_thermal_protocol_field_specific_release',reason='Internal-mode partition and matched Fc correction integrated and independently qualified; pending physical fields remain blank.',registered_additional_candidates=s['additional_candidate_count'],original39_with_native_potentials=s['original39_native_potential_molecules'],complete_potential_pairs=s['paired_complete_potential_count'],primary_candidate_snapshot=feed['snapshot_path'],arithmetic_checks=s['arithmetic_checks'],arithmetic_discrepancies=0,mode_count_audit='results/validation/rigid_mode_thermal_2026-10-18/mode_count_audit.json',independent_projected_review=str(thermal_path.relative_to(ROOT)),scientific_validation_complete=False,updated=time.time())
atomic_write_json(ROOT/'results/candidates/QUALITY_STATUS.json',quality)
print(json.dumps(dict(candidate_pairs=s['paired_complete_potential_count'],global_status=global_review['status'],slides_status=slides['status'])))
