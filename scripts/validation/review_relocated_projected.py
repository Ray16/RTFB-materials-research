#!/usr/bin/env python
"""Independently compare relocated partitions, immutable inputs and file traces."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import time


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def review(checkout,manifest,trace,forbidden,out):
    checkout=checkout.resolve();inputs=json.loads(manifest.read_text())
    for relative,expected in inputs['files_sha256'].items():
        target=(checkout/relative).resolve()
        if not target.is_relative_to(checkout) or sha(target)!=expected:
            raise ValueError('Staged input changed: '+relative)
    run=json.loads((checkout/'partition_run/run.json').read_text())
    binding=json.loads((checkout/'partition_run/runtime_binding.json').read_text())
    if run.get('checked')!=2 or run.get('skipped') or run.get('failed') or len(run['results'])!=2:
        raise ValueError('Two complete projected controls required')
    if any(r.get('kind')!='computed' for r in run['results']):
        raise ValueError('Must recompute partitions, not reuse caches')
    rows=[]
    for entry in inputs['states']:
        source=Path(entry['reference_partition'])
        if sha(source)!=entry['reference_sha256']:raise ValueError('Reference partition changed')
        reference=json.loads(source.read_text());state=checkout/'calcs/dft'/entry['id']/entry['state']
        paths=list((state/'projected_thermal').glob('*.json'))
        if len(paths)!=1:raise ValueError('One recomputed partition required')
        result=json.loads(paths[0].read_text());raw=next((state/'native_thermal').glob('*.json'))
        if result['identity']['raw_native_sha256']!=sha(raw):raise ValueError('Native bytes do not match partition')
        if result['identity']['geometry_sha256']!=sha(state/'opt.xyz'):raise ValueError('Geometry differs')
        differences={k:abs(result['nuclear_by_cutoff_eV'][k]-reference['nuclear_by_cutoff_eV'][k]) for k in ('25','50','100')}
        frequencies=result['frequencies_cm'];reference_frequencies=reference['frequencies_cm']
        if not frequencies or not all(math.isfinite(v) for v in frequencies+reference_frequencies):
            raise ValueError('Finite nonempty mode lists required')
        mode_difference=max(abs(a-b) for a,b in zip(frequencies,reference_frequencies))
        if (any(not math.isfinite(v) or v>1e-8 for v in differences.values()) or mode_difference>1e-6
                or len(result['frequencies_cm'])!=len(reference['frequencies_cm'])
                or result['symmetry_number']!=reference['symmetry_number']
                or result['implementation_certificate_sha256']!=reference['implementation_certificate_sha256']
                or abs(result['g_spin_partition_eV']-reference['g_spin_partition_eV'])>1e-12):
            raise ValueError('Relocated scientific quantities disagree')
        worker=result['kernel_execution']
        if (worker['library_sha256']!=binding['runtime']['library_sha256']
                or set(worker['thread_limits'].values())!={'1'} or len(worker['cpu_affinity'])!=1):
            raise ValueError('Kernel identity or CPU bounds differ')
        rows.append(dict(id=entry['id'],maximum_frequency_difference_cm=mode_difference,
                         nuclear_difference_eV_by_cutoff=differences,partition_sha256=sha(paths[0]),
                         native_record_sha256=sha(raw),symmetry_number=result['symmetry_number']))
    for control in binding['relocation']['controls'].values():
        local=Path(control['relocated_path']).resolve()
        if not local.is_relative_to(checkout) or sha(local)!=control['sha256']:
            raise ValueError('Relocated certificate backing bytes differ')
    text=trace.read_text();lines=text.splitlines();prefix=re.escape(str(forbidden.resolve()))
    syscall=re.compile(r'\b(open|openat|openat2|stat|lstat|statx|newfstatat|access|faccessat|faccessat2|readlink|readlinkat)\([^\n]*?"('+prefix+r'(?:/[^"\n]*)?)"')
    metadata=[]
    for line in lines:
        match=syscall.search(line)
        if not match:continue
        if match[2]==str(forbidden.resolve()) and match[1] in ('stat','lstat','newfstatat'):
            metadata.append(line)
        else:raise ValueError('Original checkout file access in traced process tree: '+line)
    xtb=binding['runtime']['binary'];cli=sum(('execve("'+xtb+'"') in line for line in lines)
    library=binding['runtime']['library'];loads=sum('openat(' in line and ('"'+library+'"') in line for line in lines)
    if cli!=6 or loads<1:raise ValueError('Missing traced CLI or explicit-library use')
    negative=json.loads((checkout/'adverse_inputs.json').read_text())
    if len(negative['checks'])!=10 or not all(c['rejected'] for c in negative['checks']):
        raise ValueError('Physical negative controls incomplete')
    record=dict(status='relocated_projected_partitions_and_physical_guards_independently_checked',checked=time.time(),
                controls=rows,immutable_staged_files=len(inputs['files_sha256']),negative_controls=10,
                traced_file_syscall_lines=len(lines),original_checkout_file_accesses=0,
                original_root_directory_metadata_calls=len(metadata),xtb_CLI_executions=cli,qualified_library_open_calls=loads,
                trace_sha256=sha(trace),runtime=binding['runtime'],relocation=binding['relocation'],
                source_sha256={str(p):sha(p) for p in (manifest,trace,checkout/'partition_run/run.json',
                               checkout/'partition_run/runtime_binding.json',checkout/'adverse_inputs.json',Path(__file__))},
                quantum_calculations_run=False,native_Hessian_execution_requalified=False,
                general_platform_portability_inferred=False)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('x') as stream:json.dump(record,stream,indent=2);stream.write('\n')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('checkout','manifest','trace','forbidden','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();record=review(args.checkout,args.manifest,args.trace,args.forbidden,args.out)
    print(json.dumps(dict(status=record['status'],controls=len(record['controls']),negative_controls=record['negative_controls'])))
