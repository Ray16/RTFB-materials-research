#!/usr/bin/env python
"""Update the maintained current collection and existing historical comparison figures."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import subprocess
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(Path(__file__).resolve().parent))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--feed',type=Path,required=True)
    ap.add_argument('--status',type=Path,required=True);a=ap.parse_args()
    import matplotlib
    matplotlib.use('Agg')
    import plot_leads_vs_existing as leads
    import plot_discovery as discovery
    leads.main(a.feed)
    discovery.main(a.feed)
    import update_comparison_collection
    update_comparison_collection.main(['--feed',str(a.feed)])
    payload=json.loads(a.feed.read_text())
    files=[ROOT/'results/figures/candidates'/name for name in [
        'leads_vs_existing_E1_disp.png','leads_vs_existing_lambda_cap.png','discovery_capacity_vs_sa.png']]
    status=dict(status='existing_plots_updated',snapshot=payload['snapshot_path'],created=time.time(),
                source_fields='admitted native potentials/dispersion; QC-clean gas diagnostic; exact-graph capacity/SA',
                previous_points_preserved=True,historical_thermal_reference_labelled=True,
                figures={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    subprocess.run([sys.executable,str(ROOT/'scripts/pipeline/publish_candidate_decisions.py'),'--feed',str(a.feed),'--comparison-already-updated'],check=True,cwd=ROOT)
    a.status.parent.mkdir(parents=True,exist_ok=True);a.status.write_text(json.dumps(status,indent=2)+'\n')


if __name__=='__main__':main()
