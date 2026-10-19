"""Maintain the recorded enumeration/selection diagram without implying all were computed."""
import hashlib,json
from collections import Counter
from pathlib import Path


def selection_overview(payload,out,root):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
    source=Path(root)/'results/discovery/systematic_families_2026-10-18/summary.json';pool=json.loads(source.read_text())
    rows=payload['additional_candidates'];matched=[r for r in rows if r['id'].startswith('search_')]
    others=[r for r in rows if not r['id'].startswith('search_')]
    counts=dict(targeted=sum(r['id'].startswith('hyp_') and not r['id'].startswith('hyp_cata_') for r in others),cata=sum(r['id'].startswith('hyp_cata_') for r in others),diaryl=sum(r['id'].startswith('study_aryl_') for r in others),long_NDI=sum(r['id']=='exp_ndi_hexyl_pentyloxy' for r in others))
    assert sum(counts.values())==len(others)
    families=Counter(r['search_core'] for r in matched)
    plt.rcdefaults();plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,ax=plt.subplots(figsize=(7.2,4.4));fig.subplots_adjust(left=.02,right=.98,top=.96,bottom=.07);ax.set(xlim=(0,100),ylim=(0,100));ax.axis('off')
    def box(x,y,w,h,title,body,accent=False):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.6,rounding_size=1.2',linewidth=.8,edgecolor='#94aab4',facecolor='#e9f4f4' if accent else '#f3f5f7'))
        ax.text(x+w/2,y+h-4,title,ha='center',va='top',fontsize=10,fontweight='bold',color='#15304a')
        ax.text(x+w/2,y+h/2-2,body,ha='center',va='center',fontsize=8,linespacing=1.55,color='#253847')
    def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=12,lw=1.1,color='#0e7c86'))
    box(2,73,44,24,f'{pool["exact_grafts"]} enumerated structures','459 viologen · 243 NDI · 135 PMDI\nSide chain × linker × ring substitution')
    box(54,73,44,24,'Literature + targeted hypotheses',f'{counts["targeted"]} targeted variants + {counts["cata"]} Cata cores\n{counts["diaryl"]} diaryl viologens + {counts["long_NDI"]} long-chain NDI')
    box(2,36,44,25,f'{len(matched)} matched comparisons','6 per family: chain length, F, CN,\nOMe, Me and linker changes')
    box(54,36,44,25,f'{len(others)} complementary designs','Core topology, attachment patterns\nand side-chain chemistry')
    arrow((24,72),(24,63));arrow((76,72),(76,63));arrow((24,35),(43,24));arrow((76,35),(57,24))
    box(19,1,62,23,f'{len(rows)} additional candidates',f'{len(payload["original39"])} original + {len(rows)} additional = {len(payload["original39"])+len(rows)} designs\nSelected calculations → field-specific QC → comparison',True)
    fig.text(.5,.015,'Enumeration and selection counts; computational completion is tracked separately.',ha='center',fontsize=7,color='#50616e')
    for suffix in ('png','pdf','svg'):fig.savefig(out/f'candidate_selection_overview.{suffix}',dpi=600)
    plt.close(fig)
    report=dict(status='recorded_family_search_and_selection',enumerated=pool['exact_grafts'],by_family=pool['by_family'],matched=len(matched),matched_by_family=dict(families),complementary=counts,additional=len(rows),original=len(payload['original39']),snapshot=payload['snapshot_path'],current_synthesis_route='post-functionalization of a preformed polymer',earlier_vinyl_branch='historical exploratory comparison; separate structures',source_sha256={str(source):hashlib.sha256(source.read_bytes()).hexdigest()})
    (out/'candidate_selection_overview.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
