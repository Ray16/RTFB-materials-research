"""Current engine policy for the reservation-aware shard launcher."""
from pathlib import Path

from redox.core.common import ROOT
from redox.backends.uma_model_identity import digest
from redox.execution import launch


def commands(module, *args, **kwargs):
    jobs = launch.commands(module, *args, **kwargs)
    for job in jobs:
        index = job.index("redox.backends." + module, job.index("--") + 1)
        job[index:index + 1] = ["redox.backends.uma_current", "--stage", module]
    return jobs


def source_metadata():
    paths = (Path(__file__).resolve(), ROOT / "src/redox/backends/uma_current.py",
             ROOT / "config/uma_execution.json")
    return {"current_entrypoint_source_sha256": {str(p): digest(p) for p in paths},
            "UMA_backend_selection": "per reserved GPU at run time: stock by default; Nitro only on a "
                                     "qualified cc 8.x-12.x card (redox.backends.uma_current.select_backend)"}


def main(argv=None):
    return launch.main(argv, command_builder=commands, report_metadata=source_metadata)


if __name__ == "__main__":
    raise SystemExit(main())
