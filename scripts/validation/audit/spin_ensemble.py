"""Independent replay of an explicit equilibrium spin selection.

No production spin-ensemble, thermal producer or spin-report helper is imported.
The caller supplies its independent single-surface/native-partition readers.
"""
import hashlib
import json
import math
from pathlib import Path


def read_ensemble(state_dir, policy, *, root, manifest, native_protocol,
                  surface_reader, partition_reader, spin_decision, strict=False):
    try:
        return _replay(state_dir, policy, root, manifest, native_protocol,
                       surface_reader, partition_reader, spin_decision)
    except (OSError, ValueError, KeyError, TypeError, AttributeError, StopIteration):
        if strict:
            raise
        return None


def _replay(state_dir, policy, root, manifest, native_protocol,
            surface_reader, partition_reader, spin_decision):
    from redox.core.protocol import ACTIVE_SP, load_record, record_path
    from redox.core.common import smiles_to_xyz_map
    sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
    root, directory=Path(root).resolve(), Path(state_dir).resolve()
    gid, state=directory.parent.name, directory.name
    if directory != (root/'calcs/dft'/gid/state).resolve():
        raise ValueError('Canonical ensemble identity differs')
    if (type(policy.get('schema')) is not int or policy['schema']!=1
            or policy.get('model')!='two-spin equilibrium; component spin degeneracy already included'
            or policy.get('temperature_K')!=298.15
            or policy.get('single_surface_kinetics_validated') is not False):
        raise ValueError('Ensemble model or kinetic claim differs')
    rows=[r for r in manifest if r['id']==gid and r['state']==state]
    if len(rows)!=1:raise ValueError('Unique declared state required')
    row=rows[0];raw=json.loads((directory/'result.json').read_text())
    q,m=raw['charge'],raw['mult']
    if type(q) is not int or type(m) is not int or q!=int(row['charge']) or not spin_decision(gid,state,raw):
        raise ValueError('Selected charge or spin differs from declaration')
    alt={1:3,3:1,2:4,4:2}[m]
    helper=(root/'calcs/spincheck'/gid/f'{state}_m{alt}').resolve()
    if policy.get('helper')!=str(helper):raise ValueError('Alternative helper differs')
    identity=dict(id=gid,state=state,charge=q,mult=m,geometry_sha256=sha(directory/'opt.xyz'),
                  smiles_sha256=hashlib.sha256(row['smiles'].encode()).hexdigest(),
                  sp_record_sha256=sha(record_path(directory,q,m)),protocol=ACTIVE_SP)
    if policy.get('chosen_identity')!=identity:raise ValueError('Selected source changed')
    pieces=[];inputs={};partition_protocol=None;seed_identity=None
    for role,folder,mult in [('chosen',directory,m),('alternative',helper,alt)]:
        data=json.loads((folder/'result.json').read_text())
        if (type(data.get('charge')) is not int or type(data.get('mult')) is not int
                or data['charge']!=q or data['mult']!=mult or data.get('optimized') is not True):
            raise ValueError('Component charge, spin or optimization differs')
        expected=dict(opt_xc='r2scan',opt_basis='def2-svpd' if q<0 else 'def2-svp',
                      opt_disp='d4',solvent='acetonitrile')
        if any(str(data.get(k)).lower()!=v for k,v in expected.items()):
            raise ValueError('Component geometry protocol differs')
        if role=='alternative':
            seed_identity=data.get('spin_source_identity')
            if (not isinstance(seed_identity,dict) or any(seed_identity.get(k)!=identity[k]
                    for k in ('id','state','charge','mult','smiles_sha256','protocol'))):
                raise ValueError('Alternative belongs to another chemical surface')
        if smiles_to_xyz_map(row['smiles'],folder/'opt.xyz') is None:
            raise ValueError('Component exact graph differs')
        tp=folder/'thermal_policy.json';thermal_selection=json.loads(tp.read_text())
        if (thermal_selection.get('thermal_definition')!='native_solution_dft_qrrho_spin'
                or thermal_selection.get('native_protocol')!=native_protocol):
            raise ValueError('Component native model is not selected')
        correction,ni,definition=surface_reader(folder,data)
        partition=partition_reader(folder,data)
        energy=load_record(folder,q,mult)
        if (correction is None or ni!=0 or definition!='native_solution_dft_qrrho_spin'
                or partition is None or energy is None):
            raise ValueError('A component is not independently admitted')
        contamination=energy.get('smd_spin_contam')
        if type(contamination) not in (float,int) or not math.isfinite(contamination) or abs(contamination)>.1:
            raise ValueError('Component spin contamination differs')
        protocol=partition['protocol']
        if protocol.get('temperature_K')!=298.15 or (partition_protocol is not None and protocol!=partition_protocol):
            raise ValueError('Component partition models differ')
        partition_protocol=protocol
        ident=dict(geometry_sha256=sha(folder/'opt.xyz'),charge=q,mult=mult,protocol=native_protocol)
        address=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16]
        native=folder/'native_thermal'/f'{address}.json';post=Path(partition['path'])
        required=[folder/'result.json',folder/'opt.xyz',tp,record_path(folder,q,mult),
                  native,native.with_suffix('')/'hessian_cartesian.npy',post]
        saved=json.loads(post.read_text())
        for name,h in saved['artifacts_sha256'].items():
            original=Path(name);local=post.with_suffix('')/original.name
            path=local if local.is_file() else original
            if sha(path)!=h:raise ValueError('Partition artifact changed')
            required.append(path)
        inputs.update({str(p):sha(p) for p in required})
        pieces.append(dict(role=role,directory=str(folder),charge=q,mult=mult,
                           electronic_eV=energy['e_smd_eV'],thermal_eV=correction,
                           free_energy_eV=energy['e_smd_eV']+correction,
                           geometry_sha256=sha(folder/'opt.xyz'),energy_path=str(record_path(folder,q,mult)),
                           native_path=str(native),partition_path=str(post)))
    if pieces[1]['electronic_eV']<pieces[0]['electronic_eV']-1e-4:
        raise ValueError('Lower electronic surface needs canonical promotion')
    saved_components=policy.get('components')
    if not isinstance(saved_components,list) or len(saved_components)!=2:
        raise ValueError('Two bound components required')
    for stored,actual in zip(saved_components,pieces):
        if set(stored)!=set(actual):raise ValueError('Component fields differ')
        for key,value in actual.items():
            if key in ('thermal_eV','free_energy_eV'):
                number=stored[key]
                if type(number) not in (int,float) or not math.isfinite(number) or abs(number-value)>1e-8:
                    raise ValueError('Independent component free energy differs')
            elif stored[key]!=value:raise ValueError('Component identity differs')
    if (policy.get('alternative_seed_identity')!=seed_identity
            or policy.get('native_partition_protocol')!=partition_protocol
            or policy.get('evidence_sha256')!=inputs):
        raise ValueError('Bound component evidence differs')
    # Factor the partition by the selected component rather than the producer's
    # minimum-energy form. Degeneracy was included by each independent native reader.
    kt=8.617333262145e-5*298.15
    gap=pieces[1]['free_energy_eV']-pieces[0]['free_energy_eV'];log_weight=-gap/kt
    log_partition=max(0.,log_weight)+math.log1p(math.exp(-abs(log_weight)))
    effective=pieces[0]['free_energy_eV']-kt*log_partition
    population=math.exp(log_weight-log_partition)
    expected=dict(free_energy_eV=effective,correction_eV=effective-pieces[0]['free_energy_eV'],
                  alternative_population=population,gap_eV=gap)
    for key,value in expected.items():
        saved=policy.get('equilibrium',{}).get(key)
        if type(saved) not in (int,float) or not math.isfinite(saved) or abs(saved-value)>1e-8:
            raise ValueError('Stored ensemble arithmetic differs')
    if any(sha(p)!=h for p,h in inputs.items()):raise ValueError('Ensemble inputs changed during replay')
    return dict(expected,components=pieces,effective_thermal_eV=effective-pieces[0]['electronic_eV'],
                selected_surface_G_eV=pieces[0]['free_energy_eV'],helper=str(helper),
                evidence_sha256=inputs,single_surface_kinetics_validated=False)
