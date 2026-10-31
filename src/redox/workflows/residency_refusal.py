"""Recover a proven foreign-residency refusal before the scientific guard starts.

This never clears GPU reservations or interrupts a peer. A namespace whose
host/scientific guard or raw outputs exist is outside this recovery path.
"""
import ast
import json
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources
from redox.workflows.submission_refusal import producer_inventory

MESSAGE='Foreign GPU residency found before scientific program start'


def refusal_branch(source,logical_path,log):
    tree=ast.parse(source.read_text());main=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    calls=[n.lineno for n in ast.walk(main) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='run_path']
    lines=log.splitlines()
    if not calls or not lines or lines[0]!='Traceback (most recent call last):' or lines[-1]!='RuntimeError: '+MESSAGE:return None
    for statement in main.body:
        if not (isinstance(statement,ast.If) and isinstance(statement.test,ast.Call)
                and isinstance(statement.test.func,ast.Name) and statement.test.func.id=='check'):continue
        for node in statement.body:
            if (isinstance(node,ast.Raise) and isinstance(node.exc,ast.Call) and isinstance(node.exc.func,ast.Name)
                    and node.exc.func.id=='RuntimeError' and len(node.exc.args)==1
                    and isinstance(node.exc.args[0],ast.Constant) and node.exc.args[0].value==MESSAGE
                    and node.lineno<min(calls) and f'  File "{logical_path}", line {node.lineno}, in main' in lines):
                return dict(message=MESSAGE,raise_line=node.lineno,first_run_path_line=min(calls))
    return None


def release_residency_refusal(spec,output,alive,hosts):
    case=Path(spec['attempt']);claim=case/'native_dispatch_claim';guard=case/'residency_guard'
    needed=[claim/'owner.json',claim/'launch.json',guard/'owner.json',guard/'state.json']
    if not all(p.is_file() for p in needed):return False
    owner,launch,guard_owner,guard_state=[json.loads(p.read_text()) for p in needed]
    if owner.get('queue')!=str(output) or launch.get('status')!='controller_spawned_scientific_start_unverified':return False
    namespace=Path(launch['namespace']);digest=launch['namespace_sha256'];args=launch['argv']
    if sha256_file(namespace)!=digest or args.count('--out')!=1:return False
    placement=Path(args[args.index('--out')+1]);closed=placement/'completion.json';log=placement/'command.log'
    if not placement.resolve().is_relative_to((output/'placements').resolve()) or not closed.is_file() or not log.is_file():return False
    terminal=json.loads(closed.read_text());logical=Path(spec['prior_root'])/'scripts/fleet/exclusive_python_worker.py'
    source=Path(spec['archive'])/'scripts/fleet/exclusive_python_worker.py'
    branch=refusal_branch(source,logical,log.read_text())
    command=terminal.get('argv',[])
    archive_bind=any(command[i:i+3]==['--ro-bind',spec['archive'],spec['prior_root']] for i in range(len(command)-2))
    if (branch is None or not archive_bind or terminal.get('returncode') in (None,0) or terminal.get('namespace')!=str(namespace)
            or terminal.get('namespace_sha256')!=digest or str(logical) not in terminal.get('argv',[])
            or guard_owner['host'].split('.')[0]!=terminal['host'].split('.')[0]
            or guard_owner['GPU']!=terminal['GPU'] or guard_state['GPU']!=terminal['GPU']
            or not guard_state.get('foreign_resident_PIDs') or guard_state.get('own_residency_observed') is not False):return False
    absent=[case/name for name in ('host_guard','tensor_runtime','experimental_native','effective_execution.json','independent_native_review.json')]
    absent+=[NamespaceSources(spec).path(spec['original_raw_root'])]
    if any(p.exists() for p in absent):return False
    owners=[launch['owner'],guard_owner,dict(host=terminal['host'],pid=terminal['pid'],process_start=terminal['process_start'])]
    proof_paths=[*needed,closed,log,source,namespace,case/'plan.json',case/'input.xyz']
    before={str(p):sha256_file(p) for p in proof_paths}
    try:
        if any(alive(o) is not False for o in owners):return False
        inventory=producer_inventory(spec,namespace,hosts)
        if set(inventory)!=set(hosts) or any(inventory.values()):return False
        if any(p.exists() for p in absent) or any(alive(o) is not False for o in owners):return False
        if any(sha256_file(p)!=h for p,h in before.items()):return False
    except (OSError,ValueError,subprocess.SubprocessError):return False
    history=output/'residency_refusal_history'/spec['key']/str(time.time_ns());history.mkdir(parents=True)
    # The queue holds its exclusive lock throughout. No matching producer can
    # own this namespace; preserve the complete guard and claim as history.
    shutil.move(str(guard),str(history/'residency_guard'));shutil.move(str(claim),str(history/'claim'))
    copied=history/'exclusive_worker_source.py';shutil.copyfile(source,copied)
    preserved={}
    for path,digest in before.items():
        path=Path(path)
        if path.is_relative_to(guard):path=history/'residency_guard'/path.relative_to(guard)
        elif path.is_relative_to(claim):path=history/'claim'/path.relative_to(claim)
        preserved[str(path)]=digest
    if any(sha256_file(p)!=h for p,h in preserved.items()):raise ValueError('Refusal evidence changed during preservation')
    atomic_write_json(history/'resolution.json',dict(status='pre_science_residency_refusal_released_without_scientific_attempt',
        namespace=str(namespace),namespace_sha256=digest,branch=branch,stopped_owners=owners,
        producer_inventory=inventory,placement_preserved=str(placement),scientific_outputs_absent=True,
        scientific_attempt_consumed=False,foreign_processes_signalled=False,source_sha256=preserved,
        original_source_sha256=before,
        immutable_guard_source=str(copied),immutable_guard_source_sha256=sha256_file(copied),checked=time.time()))
    return True
