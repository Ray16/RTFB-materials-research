"""Inventory the frozen publication component records and declared dependencies."""
import json
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.validation.source_inventory import inventory
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

with allocation_guard():
    host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
ROOT=Path(__file__).resolve().parents[2]
FOLDER=Path(__file__).resolve().parent
OUT=ROOT.parent/'redox_release_candidates/2026-11-03_publication_raw_inventory.json'
if OUT.exists():raise FileExistsError('Preserve earlier inventory')
seed_path=FOLDER/'raw_seed_manifest.json';seeds=json.loads(seed_path.read_text())
known_path=ROOT/'validation/counterion_release_inventory_2026-11-03/inventory.json'
known=json.loads(known_path.read_text())['files']
extra_path=ROOT/'validation/counterion_release_inventory_2026-11-03/grafting_metadata_resolution.json'
known.extend(json.loads(extra_path.read_text())['resolved'])
report=inventory(ROOT,ROOT.parent/'archive/0_redox',seeds['seeds_sha256'],known)
report['source_sha256']={str(p):sha256_file(p) for p in
    (seed_path,known_path,extra_path,Path(__file__),Path(inventory.__code__.co_filename))}
OUT.write_text(json.dumps(report,indent=2)+'\n')
summary=dict(status=report['status'],inventory=str(OUT),inventory_sha256=sha256_file(OUT),
    verified_files=len(report['files']),unresolved_bindings=len(report['unresolved']),
    unique_content_bytes=report['unique_content_bytes'],seed_groups=seeds['groups'],
    implicit_dependencies_covered=False,quantum_validation_repeated=False,public_deposition_complete=False)
(FOLDER/'raw_inventory_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='seed_groups'}))
