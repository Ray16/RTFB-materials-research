"""Matched-cutoff sensitivity of the new equilibrium result; 50 cm-1 stays primary."""
from pathlib import Path
import hashlib,math,os,time
from redox.core.common import ROOT
from redox.core.protocol import load_record,record_path
from redox.core.io import atomic_write_json
from redox.qm.native_thermal import load_native_thermal,cache_path
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
base=Path(__file__).resolve().parent;gid='search_pmdi_1864c2c74267'
specs=[('red1',ROOT/'calcs/dft'/gid/'red1',-1,2),('singlet',ROOT/'calcs/dft'/gid/'red2',-2,1),('triplet',ROOT/'calcs/spincheck'/gid/'red2_m3',-2,3),('Fc',ROOT/'calcs/dft/ferrocene/neu',0,1),('Fcox',ROOT/'calcs/dft/ferrocene/ox',1,2)]
data={};evidence={};sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for key,path,q,m in specs:
 sp=load_record(path,q,m);th=load_native_thermal(path,q,m);assert sp and th
 data[key]=(sp['e_smd_eV'],th)
 for p in [path/'opt.xyz',record_path(path,q,m),cache_path(path,q,m)[0],Path(th['thermal_postprocessing_path'])]:evidence[str(p)]=sha(p)
kt=8.617333262145e-5*298.15;rows=[]
for cutoff in [25,50,100]:
 G={k:e+t['nuclear_by_cutoff_eV'][str(cutoff)]+t['g_spin_partition_eV'] for k,(e,t) in data.items()}
 gap=G['triplet']-G['singlet'];logweight=-gap/kt;logsum=max(0.,logweight)+math.log1p(math.exp(-abs(logweight)))
 mixed=G['singlet']-kt*logsum;fc=G['Fcox']-G['Fc']
 rows.append(dict(cutoff_cm=cutoff,triplet_population=math.exp(logweight-logsum),G_mixture_eV=mixed,Fc_absolute_V=fc,E2_vs_Fc_V=G['red1']-mixed-fc))
for p,h in evidence.items():assert sha(p)==h
result=dict(status='matched_cutoff_equilibrium_sensitivity',primary_cutoff_cm=50,rows=rows,E2_range_V=max(r['E2_vs_Fc_V'] for r in rows)-min(r['E2_vs_Fc_V'] for r in rows),source_sha256=sha(Path(__file__)),evidence_sha256=evidence,production_protocol_changed=False,experimental_accuracy_validated=False,finished=time.time())
atomic_write_json(base/'cutoff_sensitivity.json',result);print(rows);print('E2 range',result['E2_range_V'])
