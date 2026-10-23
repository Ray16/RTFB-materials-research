#!/usr/bin/env python
"""Exact-candidate adjacent-site feasibility searches, linked to enumeration.

Existing capped-molecule redox results are retained. These neutral resting-state
polymer searches do not constitute a new full-redox campaign or synthesis pass.
The common empty backbone is identified explicitly instead of computed again.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from rdkit import Chem
from rdkit.Chem import AllChem

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.build.polymer_cases import build_case
from redox.build.structures import to_xyz
from redox.build.geometry import geometry_stereo_status
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def prepare(out, identifier):
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    feed_path = ROOT / "results/candidates/current.json"
    feed = json.loads(feed_path.read_text())
    candidate = next(row for row in feed["original39"] + feed["additional_candidates"] if row["id"] == identifier)
    evidence_path = ROOT / "results/grafting/route_evidence.json"
    evidence = json.loads(evidence_path.read_text())
    route_module = api.load_module("candidate_occupancy_route_gate", ROOT / "scripts/pipeline/screen_grafting_routes.py")
    assessment = json.loads((ROOT / "results/grafting/selection_assessment.json").read_text())
    route = route_module.screen(candidate, evidence["candidates"], assessment)
    if not route["proceed_to_feasibility_model"]:
        raise ValueError("Prospective route gate excludes feasibility modeling")
    if evidence["candidates"].get(identifier, {}).get("exact_precursor_access_reviewed") is not True:
        raise ValueError("This focused end-to-end campaign requires a reviewed exact precursor")
    key = Chem.MolToInchiKey(Chem.MolFromSmiles(candidate["exact_grafted_smiles"]))
    pool_path = ROOT / "results/discovery/systematic_families_2026-10-18/structures.json"
    pool = json.loads(pool_path.read_text())
    aliases = [row["group"]["id"] for row in pool if Chem.MolToInchiKey(Chem.MolFromSmiles(row["grafted_smiles"])) == key]
    if not aliases:
        raise ValueError("Selected candidate does not exactly match the enumerated pool")
    out.mkdir(parents=True, exist_ok=False)
    atomic_write_json(out / "candidate_snapshot.json", candidate)
    atomic_write_json(out / "route_snapshot.json", dict(sources=evidence["sources"], candidate=evidence["candidates"][identifier], gate=route))
    fragment = candidate["frag"]
    empty = build_case(fragment, 5, [2, 3, 4], [], ["CW"] * 4)
    shared_path = ROOT / "results/validation/conformer_sampling_2026-10-23/reference_repeats/all_same/case.json"
    shared = json.loads(shared_path.read_text())
    if empty["smiles"] != shared["smiles"] or empty["backbone_stereo_fingerprint"] != shared["backbone_stereo_fingerprint"]:
        raise ValueError("Proposed common reference is not the same constitutional/stereo backbone")
    tasks = []
    for occupied in ([2], [3], [2, 3]):
        case = build_case(fragment, 5, [2, 3, 4], occupied, ["CW"] * 4)
        if case["charge"] != 0:
            raise ValueError("Charged motifs require a separately qualified ionic model")
        label = identifier + "__occ_" + "_".join(map(str, occupied))
        directory = out / label
        directory.mkdir()
        mol = Chem.AddHs(Chem.MolFromSmiles(case["smiles"]))
        parameters = AllChem.ETKDGv3()
        parameters.randomSeed = 20261023
        parameters.numThreads = 1
        conformers = list(AllChem.EmbedMultipleConfs(mol, numConfs=32, params=parameters))
        optimized = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=2000)
        converged = [(energy, cid) for cid, (status, energy) in zip(conformers, optimized) if status == 0]
        if not converged:
            raise ValueError("No converged MMFF starting geometry for " + label)
        energy, cid = min(converged)
        one = Chem.Mol(mol)
        one.RemoveAllConformers()
        one.AddConformer(mol.GetConformer(cid), assignId=True)
        xyz = directory / "input.xyz"
        xyz.write_text(to_xyz(one, "MMFF warm start only; exact-candidate polymer feasibility search follows"))
        if hydrogen_topology_status(case["smiles"], xyz) != "ok" or geometry_stereo_status(case["smiles"], xyz) not in ("ok", "not_applicable"):
            raise ValueError("Candidate occupancy input failed graph/stereo QC")
        case.update(key=label, mult=1, motif_id=identifier, fragment=fragment, source_capped_product_key=key,
            initial_embedding_seed=20261023, embedded=len(conformers), converged_MMFF=len(converged),
            initial_MMFF_energy_kcal_mol=energy, role="Exact-candidate occupancy feasibility model; no synthesis admission")
        atomic_write_json(directory / "case.json", case)
        tasks.append(dict(key=label, input=str(xyz), case=str(directory / "case.json"),
            raw_out=str(ROOT / "calcs/polymer_models" / out.name / label)))
    sources = [Path(__file__).resolve(), pool_path, out / "candidate_snapshot.json", out / "route_snapshot.json", shared_path]
    atomic_write_json(out / "plan.json", dict(candidate_id=identifier, exact_capped_product_key=key, enumerated_aliases=aliases,
        tasks=tasks, common_empty_reference=dict(case_path=str(shared_path), graph_smiles=empty["smiles"],
            broader_search="calcs/polymer_models/2026-10-23_conformer_revision/all_same__occ_none__parallel_repeat",
            baseline_native_campaign="results/polymer_models/campaigns/2026-10-23_baseline_two_pose_ensemble_v2",
            scope="Same empty graph only; no transfer of occupied phthalimide energies or uncomputed ensemble G"),
        source_sha256={str(p): api.digest(p) for p in sources}, created=time.time(),
        source_feed_as_of=feed["as_of_Chicago"], precursor_access_reviewed=True, prospective_route_status=route["route_status"],
        full_redox_campaign_submitted=False, synthesis_admission=False,
        scope="Three adjacent-site occupancy graphs for one exact enumerated PMDI design; common empty component identified separately"))


def launch(out, host):
    preflight = json.loads((out / (host + "_preflight.json")).read_text())
    if preflight.get("ok") is not True or preflight["load"] + 42 > preflight["online_cpus"]:
        raise ValueError("Node preflight cannot cover three 14-core source searches")
    plan = json.loads((out / "plan.json").read_text())
    api.verify_sources(plan["source_sha256"])
    remote = '''import os,json,subprocess,time,sys
from pathlib import Path
c=json.load(sys.stdin);rows=[]
for task in c['tasks']:
 if Path(task['raw_out']).exists():raise ValueError('Existing search must be observed, not restarted')
 env=dict(os.environ,PYTHONPATH=c['root']+'/src',OMP_NUM_THREADS='14',OMP_THREAD_LIMIT='14',MKL_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',XTB_THREADS='1')
 argv=[c['python'],c['root']+'/scripts/pipeline/sample_polymer_conformers.py','--input',task['input'],'--case',task['case'],'--out',task['raw_out'],'--threads','14','--max-hours','48']
 with open(c['out']+'/'+task['key']+'_launch.log','x') as log:p=subprocess.Popen(argv,cwd=c['root'],env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 rows.append(dict(task,pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv))
 time.sleep(.3)
print(json.dumps(dict(host=os.uname().nodename,launches=rows)))'''
    payload = dict(root=str(ROOT), python=sys.executable, out=str(out), tasks=plan["tasks"])
    result = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", host, shlex.join([sys.executable, "-c", remote])],
        input=json.dumps(payload), capture_output=True, text=True, timeout=30, check=True)
    proof = json.loads(result.stdout)
    atomic_write_json(out / "launch.json", proof)
    print(json.dumps(dict(host=host,submitted=len(proof["launches"]),actual_sampling_starts_unconfirmed=True)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "launch"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--candidate", default="design_pmdi_ethyl")
    parser.add_argument("--host", default="lambda13")
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.out.resolve(), args.candidate)
    else:
        launch(args.out.resolve(), args.host)


if __name__ == "__main__":
    main()
