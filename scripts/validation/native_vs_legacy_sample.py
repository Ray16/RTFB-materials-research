import os, sys, subprocess, tempfile, re
from pathlib import Path
import numpy as np
sys.path.insert(0, "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/src")
from ase.io import read as aseread
from redox.qm import qrrho as Q
XTB="/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb"
ARCH="/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/calcs/dft"
EV=27.211386245988; Z={"H":1,"C":6,"N":7,"O":8,"F":9,"Fe":26,"S":16}
def run_state(gid,state,q,mult):
    geom=f"{ARCH}/{gid}/{state}/opt.xyz"
    if not Path(geom).is_file(): return None
    at=aseread(geom)
    with tempfile.TemporaryDirectory() as d:
        aseread(geom).write(str(Path(d)/"mol.xyz"))
        env=dict(os.environ,OMP_NUM_THREADS="4",XTB_THREADS="4")
        p=subprocess.run([XTB,"mol.xyz","--gfn","2","--chrg",str(q),"--uhf",str(mult-1),"--hess","--acc","1.0"],
                         cwd=d,env=env,capture_output=True,text=True)
        m=re.search(r"G\(RRHO\) contrib\.\s+([-\d.]+)",p.stdout)
        lg=float(m.group(1))*EV if m else None
        hf=Path(d)/"hessian"
        if not hf.is_file() or lg is None: return None
        body=hf.read_text().split("$hessian",1)[1].split("$end",1)[0]
        vals=[float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?",body)]
        n3=3*len(at); H=np.array(vals[:n3*n3]).reshape(n3,n3)
        try: ng=Q.nuclear_thermals(H,at.positions,[Z[s] for s in at.get_chemical_symbols()],1)["nuclear_by_cutoff_eV"]["50"]
        except Exception as e: return ("err",str(e)[:60])
        return (lg,ng)
CASES={
 "quinone(NQ) d3_05stji":[("neu",0,1),("red1",-1,2),("red2",-2,1)],
 "AQ exp_aq_1oxy":[("neu",0,1),("red1",-1,2),("red2",-2,1)],
 "PMDI search_pmdi_87840c94153a":[("neu",0,1),("red1",-1,2),("red2",-2,1)],
 "NDI search_ndi_573f69a7895b":[("neu",0,1),("red1",-1,2),("red2",-2,1)],
 "viologen search_viologen_1ad45438b9e4":[("neu",0,1),("red1",-1,2),("red2",-2,1)],
}
GIDS={"quinone(NQ) d3_05stji":"d3_05stji","AQ exp_aq_1oxy":"exp_aq_1oxy",
 "PMDI search_pmdi_87840c94153a":"search_pmdi_87840c94153a","NDI search_ndi_573f69a7895b":"search_ndi_573f69a7895b",
 "viologen search_viologen_1ad45438b9e4":"search_viologen_1ad45438b9e4"}
for lab,states in CASES.items():
    gid=GIDS[lab]; print(f"\n=== {lab} ===",flush=True); G={}
    for st,q,mult in states:
        r=run_state(gid,st,q,mult)
        if r is None: print(f"  {st}: (no data)"); continue
        if r[0]=="err": print(f"  {st}: ERR {r[1]}"); continue
        G[st]=r; print(f"  {st}: legacy={r[0]:.4f} native={r[1]:.4f} diff={r[1]-r[0]:+.4f} eV",flush=True)
    order=[s for s,_,_ in states if s in G]
    for a,b in zip(order,order[1:]):
        dl=G[b][0]-G[a][0]; dn=G[b][1]-G[a][1]
        print(f"  couple {a}->{b}: basis-change potential shift = {(dn-dl)*1000:+.1f} mV",flush=True)
