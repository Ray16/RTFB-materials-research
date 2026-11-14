"""Per-state GPU fan-out of the motif DFT stage."""
import sys
import threading

import pytest

from redox.workflows import motifs


def test_state_jobs_select_one_state_each_and_keep_every_other_argument():
    argv = ["py", "-m", "redox.qm.uma_current", "--stage", "dft", "--only", "g", "--nconf", "3"]
    jobs = motifs.state_jobs(argv, "g", ["neu", "red1"])
    assert jobs == [("neu", [*argv[:6], "g:neu", "--nconf", "3"]),
                    ("red1", [*argv[:6], "g:red1", "--nconf", "3"])]
    with pytest.raises(ValueError):
        motifs.state_jobs(argv, "other", ["neu"])


def test_parallel_states_never_share_a_slot_and_a_failure_does_not_stop_the_others(tmp_path):
    lock, active, seen = threading.Lock(), set(), []

    def wrap(argv, gpu, cpus):
        with lock:
            assert gpu not in active, "two jobs on one GPU"
            active.add(gpu)
            seen.append((argv[-1], gpu, tuple(cpus)))
        return argv

    original = motifs.run_stage

    def run_stage(argv, log, env):
        try:
            return original(argv, log, env)
        finally:
            with lock:
                active.discard(next(g for s, g, _ in seen if s == argv[-1]))

    motifs.run_stage = run_stage
    try:
        code = "import sys,time; time.sleep(0.2); sys.exit(sys.argv[1] == 'red1')"
        jobs = [(s, [sys.executable, "-c", code, s]) for s in ("neu", "red1", "red2", "ox")]
        results = motifs.run_parallel(jobs, [(0, [0, 1]), (1, [2, 3])], tmp_path / "02_dft.log", None, wrap)
    finally:
        motifs.run_stage = original
    assert [state for state, _, _ in results] == ["neu", "red1", "red2", "ox"]
    assert [proc.returncode for _, _, (proc, _) in results] == [0, 1, 0, 0]
    assert {gpu for _, gpu, _ in seen} == {0, 1} and len(seen) == 4
    assert all(cpus == {0: (0, 1), 1: (2, 3)}[gpu] for _, gpu, cpus in seen)
    assert (tmp_path / "02_dft_red2.log").exists()
