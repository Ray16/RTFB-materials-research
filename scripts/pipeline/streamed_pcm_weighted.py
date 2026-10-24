"""Compose the controlled streamed curvature with tiled weighted response."""
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def compact_smd_execution(directory, budget_GiB):
    from redox.core.common import ROOT
    from redox.workflows import polymer_campaign as api
    streamed = api.load_module('retained_streamed_curvature',
        ROOT / 'scripts/pipeline/streamed_pcm_execution.py')
    weighted = api.load_module('streamed_weighted_first_response',
        ROOT / 'scripts/pipeline/smd_weighted_response.py')
    with streamed.compact_smd_execution(directory, budget_GiB=budget_GiB) as proof:
        with weighted.weighted_response(directory) as metrics:
            proof['source_sha256'].update(metrics['source_sha256'])
            proof['source_sha256'][str(Path(__file__).resolve())] = api.digest(__file__)
            proof.update(weighted_PCM_response=True, dense_weighted_surface_matrices=False)
            yield proof
