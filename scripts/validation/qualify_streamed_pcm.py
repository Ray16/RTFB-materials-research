#!/usr/bin/env python
"""Full singlet/doublet matrix replay for streamed PCM curvature execution."""
import hashlib
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def checked_controls(folder):
    evidence={}
    for mult in (1,2):
        path=Path(folder)/'controls'/str(mult)/'status.json'
        result=json.loads(path.read_text())
        if (result['status']!='full_original_compact_equivalence_passed'
                or not result.get('streamed_operator_checks_passed') or result['mult']!=mult):
            raise ValueError('Both full spin/matrix streamed-operator controls must pass')
        api.verify_sources(result['source_sha256'])
        api.verify_sources(result['evidence_sha256'])
        if any(not 0<=v<=1e-6 for v in result['differences'].values()):
            raise ValueError('The full spin-control comparisons exceed their original bounds')
        evidence[str(path)]=api.digest(path)
        evidence.update(result['source_sha256'])
    return evidence


def main():
    source=ROOT/'scripts/validation/qualify_smd_workspace.py'
    qualifier=api.load_module('streamed_PCM_full_spin_replay',source)
    before=qualifier.adapted_body
    def adapted(body):
        text=before(body)
        old='scripts/pipeline/smd_workspace_execution.py'
        if text.count(old)!=1:
            raise ValueError('The fixed full spin-control template changed')
        return text.replace(old,'scripts/pipeline/streamed_pcm_execution.py')
    qualifier.adapted_body=adapted
    try:
        qualifier.main()
    finally:
        out=Path(sys.argv[sys.argv.index('--out')+1])
        path=out/'status.json'
        if path.exists():
            result=json.loads(path.read_text())
            for p in (Path(__file__).resolve(),source,ROOT/'scripts/pipeline/streamed_pcm_execution.py'):
                result.setdefault('source_sha256',{})[str(p)]=api.digest(p)
            result.update(execution_variant='Streamed FP64 switch contractions and bounded GPU offdiagonal primitives',
                large_resource_qualified=False,native_G_production_qualified=False)
            path.write_text(json.dumps(result,indent=2)+'\n')
    stats_path=out/'compact/execution/streamed_PCM.json'
    stats=json.loads(stats_path.read_text())
    if (stats['status']!='streamed_PCM_finished_requires_full_replay_and_physical_QC'
            or min(stats[k] for k in ('switch_contractions','switch_primitive_calls','diagonal_contractions','offdiagonal_primitive_calls'))<1
            or stats['maximum_offdiagonal_primitive_bytes']>512*1024**2
            or stats['peak_GPU_block_cache_bytes']>128*1024**2):
        raise ValueError('The full replay did not exercise every bounded operator')
    api.verify_sources(stats['source_sha256'])
    result=json.loads(path.read_text())
    result['streamed_operator_checks_passed']=True
    result['source_sha256'][str(stats_path)]=api.digest(stats_path)
    path.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
