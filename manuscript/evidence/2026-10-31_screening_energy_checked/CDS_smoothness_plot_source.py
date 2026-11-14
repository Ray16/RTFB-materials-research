#!/usr/bin/env python
"""Plot measured radial CDS curvatures; no full-energy or G inference."""
import argparse
import hashlib
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--legacy',type=Path,required=True);parser.add_argument('--smooth',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    legacy=json.loads(args.legacy.read_text());smooth=json.loads(args.smooth.read_text())
    a=legacy['centered_curvature_Ha_Bohr2_by_displacement_A']
    b=smooth['states']['contact']['contact_scan']['centered_curvature_Ha_Bohr2']
    plt.rcParams.update({'font.size':9,'axes.labelsize':9,'legend.fontsize':8,'pdf.fonttype':42,
        'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    figure,axis=plt.subplots(figsize=(3.5,2.7),layout='constrained')
    plotted={}
    for name,values,color,marker in [('Legacy CDS',a,'#D55E00','s'),('Smooth CDS control',b,'#0072B2','o')]:
        points=sorted((float(k),v) for k,v in values.items());plotted[name]=points
        axis.loglog([p[0] for p in points],[p[1] for p in points],label=name,color=color,marker=marker,ms=4,lw=1.2)
    axis.set(xlabel=r'Radial displacement $h$ ($\mathrm{\AA}$)',
        ylabel=r'CDS curvature ($E_h\,a_0^{-2}$)')
    axis.legend(frameon=False,loc='center left')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    paths=[]
    for suffix in ('pdf','png'):
        p=args.out.with_suffix('.'+suffix);figure.savefig(p,dpi=300);paths.append(p)
    plt.close(figure)
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    report=dict(status='cached_CDS_component_curvature_figure_generated',data=plotted,
        source_sha256={str(p):sha(p) for p in (args.legacy,args.smooth,Path(__file__).resolve())},
        output_sha256={str(p):sha(p) for p in paths},
        scope='One radial component scan at the tested 4.2 Å carbon contact; not a full DFT Hessian or free energy')
    args.out.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],outputs=[str(p) for p in paths])))


if __name__=='__main__':main()
