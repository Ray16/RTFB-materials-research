"""Process environments shared by orchestration, without changing os.environ."""
import os
from pathlib import Path

THREAD_KEYS = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS",
               "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")


def worker_environment(root, *, threads=1, base=None, python=None):
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    env = dict(os.environ if base is None else base)
    # An inherited editable checkout must not win over the requested source tree.
    env["PYTHONPATH"] = str(Path(root) / "src")
    env.update({key: str(threads) for key in THREAD_KEYS})
    if python is not None:
        env["PATH"] = str(Path(python).parent) + os.pathsep + env.get("PATH", "")
    return env
