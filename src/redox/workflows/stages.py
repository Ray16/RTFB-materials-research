"""Dependency-aware execution of independent subprocess stages.

Only the coordinator changes status. Once a stage fails, no more work starts;
already-running independent stages finish and retain their logs.
"""
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait


def validate_graph(stages):
    names = [stage.name for stage in stages]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate stage names")
    completed = set()
    pending = list(stages)
    while pending:
        ready = [stage for stage in pending if set(stage.dependencies) <= completed]
        if not ready:
            raise ValueError("Stage graph contains a cycle or missing dependency")
        completed.update(stage.name for stage in ready)
        pending = [stage for stage in pending if stage not in ready]


def execute_graph(stages, run, *, workers=1, started=lambda stage: None,
                  finished=lambda stage, result: None):
    """Call run(stage, slot) at most once, after every prerequisite succeeds.

    run returns an integer exit code. Slot numbers let subprocess callers assign
    disjoint CPU affinity. Callbacks run only in the coordinator thread.
    """
    if type(workers) is not int or workers < 1:
        raise ValueError("workers must be a positive integer")
    validate_graph(stages)
    pending, completed, failures = list(stages), set(), []
    available = list(range(workers))
    active = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        while pending or active:
            if not failures:
                for stage in list(pending):
                    if not available:
                        break
                    if not set(stage.dependencies) <= completed:
                        continue
                    slot = available.pop(0)
                    pending.remove(stage)
                    started(stage)
                    active[pool.submit(run, stage, slot)] = (stage, slot)
            if not active:
                break
            done, _ = wait(active, return_when=FIRST_COMPLETED)
            # Consume all completions before scheduling successors or replacements.
            for future in sorted(done, key=lambda f: stages.index(active[f][0])):
                stage, slot = active.pop(future)
                available.append(slot)
                code = future.result()
                finished(stage, code)
                if code:
                    failures.append(stage.name)
                else:
                    completed.add(stage.name)
    return {"completed": sorted(completed), "failed": failures,
            "not_started": [stage.name for stage in pending]}
