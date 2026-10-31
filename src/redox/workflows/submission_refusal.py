"""Recognize resource refusals that occurred before any scientific subprocess.

Neither a dead controller nor a traceback alone authorizes reopening a task.
The exact source branch, source hash, absent outputs and fleet-wide absence of
matching producers must all agree. Unknown ownership preserves the claim.
"""
import ast
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import shlex
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources


STAGES = {
    'native': dict(claim='native_dispatch_claim', module='derived_native_namespace.py',
        messages=('GPU occupied or reserved; no native calculation started',
                  'Host RAM already promised; no native calculation started')),
    'energy': dict(claim='dispatch_claim', module='archived_energy_namespace.py',
        messages=('Device occupied or reserved; no scientific attempt started',
                  'RAM already promised to other live workers; no scientific attempt started')),
}


def refusal_branch(source, log, messages):
    """A unique literal RuntimeError must precede every Popen in launch()."""
    tree = ast.parse(source.read_text())
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'launch')
    spawns = [node.lineno for node in ast.walk(function) if isinstance(node, ast.Call)
              and isinstance(node.func, ast.Attribute) and node.func.attr == 'Popen']
    if not spawns:
        return None
    lines = log.splitlines()
    if not lines or lines[0] != 'Traceback (most recent call last):':
        return None
    for message in messages:
        raises = [node for node in ast.walk(function) if isinstance(node, ast.Raise)
                  and isinstance(node.exc, ast.Call) and isinstance(node.exc.func, ast.Name)
                  and node.exc.func.id == 'RuntimeError' and len(node.exc.args) == 1
                  and isinstance(node.exc.args[0], ast.Constant) and node.exc.args[0].value == message]
        if (len(raises) == 1 and raises[0].lineno < min(spawns)
                and lines[-1] == 'RuntimeError: '+message
                and f'  File "{source}", line {raises[0].lineno}, in launch' in lines):
            return dict(message=message, raise_line=raises[0].lineno, first_subprocess_line=min(spawns))
    return None


def producer_inventory(spec, namespace, hosts):
    needles = [str(namespace), spec['original_campaign'], spec['attempt']]
    code = ('import pathlib,json;needles='+repr(needles)+';rows=[]\n'
            'for path in pathlib.Path("/proc").glob("[0-9]*/cmdline"):\n'
            ' try:\n'
            '  args=[v.decode(errors="replace") for v in path.read_bytes().split(b"\\0") if v]\n'
            '  if any(value in needles for value in args):\n'
            '   stat=(path.parent/"stat").read_text().rsplit(")",1)[1].split()\n'
            '   if stat[0]!="Z":rows.append(dict(pid=int(path.parent.name),process_start=stat[19],argv=args))\n'
            ' except FileNotFoundError:pass\n'
            'print(json.dumps(rows))')
    def inspect(host):
        result = subprocess.run(['ssh', '-S', 'none', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host,
            shlex.join(['/usr/bin/python3', '-c', code])], capture_output=True, text=True, timeout=25, check=True)
        rows = json.loads(result.stdout)
        if not isinstance(rows, list):
            raise ValueError('Unknown producer inventory')
        return host, rows
    with ThreadPoolExecutor(max_workers=max(1, len(hosts))) as pool:
        return dict(pool.map(inspect, hosts))


def release_prelaunch_refusal(spec, output, alive, hosts, stage, source_hashes):
    policy = STAGES[stage]
    case = Path(spec['attempt'])
    claim = case/policy['claim']
    if not (claim/'owner.json').exists() or not (claim/'launch.json').exists():
        return False
    owner = json.loads((claim/'owner.json').read_text())
    launch = json.loads((claim/'launch.json').read_text())
    if stage == 'energy':
        if len(launch.get('receipts', [])) != 1:
            return False
        launch = launch['receipts'][0]
    if owner.get('queue') != str(output) or launch.get('status') != 'controller_spawned_scientific_start_unverified':
        return False
    namespace = Path(launch.get('namespace', launch.get('spec', '')))
    digest = launch.get('namespace_sha256', launch.get('spec_sha256'))
    if not namespace.is_absolute() or sha256_file(namespace) != digest:
        return False
    args = launch['argv']
    if args.count('--out') != 1:
        return False
    placement = Path(args[args.index('--out')+1])
    if not placement.resolve().is_relative_to((output/'placements').resolve()):
        return False
    log = placement.parent/'controller.log'
    source = ROOT/'src/redox/workflows'/policy['module']
    if not log.exists() or source_hashes.get(str(source)) != sha256_file(source):
        return False
    branch = refusal_branch(source, log.read_text(), policy['messages'])
    if branch is None:
        return False
    absent = [case/name for name in ('residency_guard', 'host_guard', 'tensor_runtime',
                                    'experimental_native', 'experimental_energy')]
    absent += [placement/name for name in ('launch.json', 'completion.json', 'command.log', 'launcher_source.py')]
    absent += [NamespaceSources(spec).path(spec['original_raw_root'])]
    if any(p.exists() for p in absent) or Path(spec['key']).name != spec['key']:
        return False
    try:
        if alive(launch['owner']) is not False:
            return False
        inventory = producer_inventory(spec, namespace, hosts)
        if set(inventory) != set(hosts) or any(inventory.values()):
            return False
    except (OSError, ValueError, subprocess.SubprocessError):
        return False
    # Recheck after the distributed observation before changing ownership.
    try:
        if (any(p.exists() for p in absent) or alive(launch['owner']) is not False
                or sha256_file(source) != source_hashes[str(source)] or sha256_file(namespace) != digest):
            return False
    except (OSError, ValueError, subprocess.SubprocessError):
        return False
    history = output/'prelaunch_refusal_history'/spec['key']/str(time.time_ns())
    history.mkdir(parents=True)
    shutil.move(str(claim), str(history/'claim'))
    atomic_write_json(history/'resolution.json', dict(status='prelaunch_resource_refusal_no_scientific_attempt',
        stage=stage, branch=branch, controller=launch['owner'], producer_inventory=inventory,
        source=str(source), source_sha256=source_hashes[str(source)], namespace=str(namespace),
        namespace_sha256=digest, checked=time.time(), scientific_attempt_consumed=False,
        scientific_outputs_absent=True, placement_preserved=str(placement)))
    return True
