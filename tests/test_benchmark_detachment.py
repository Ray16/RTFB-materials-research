"""Detachment repair changes only the intended remote launch, preserving quoting."""

import importlib.util
from pathlib import Path
import shlex


def test_detachment_removes_compound_subshell_and_preserves_argv():
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "benchmark_detachment", root / "scripts/pipeline/launch_native_benchmark.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    command = (
        "cd " + shlex.quote(str(root)) + " && nohup env KEY='a b' python "
        "-m redox.workflows.benchmark_native worker --campaign '/a b' "
        "> '/log a' 2>&1 < /dev/null & echo $!"
    )
    argv = ["ssh", "-n", "node", command]
    changed = module.detached_argv(argv)
    assert changed[:-1] == argv[:-1]
    assert changed[-1] == command.split(" && ", 1)[1]
    assert "&&" not in changed[-1]
    for other in (
        ["echo", command],
        ["ssh", "node", "nvidia-smi"],
        ["ssh", "node", "nohup other &"],
    ):
        assert module.detached_argv(other) == other
