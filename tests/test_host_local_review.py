import ast
import os
import shlex

from redox.workflows.host_local_review import command_for_host


def test_local_review_retains_exact_original_namespace_command():
    command = ['bwrap', '--ro-bind', '/archive', '/old', '--', 'python', 'review.py']
    assert command_for_host(command, os.uname().nodename) == command


def test_remote_review_allocates_cpu_before_entering_original_namespace():
    command = ['/shared/bwrap', '--ro-bind', '/archive', '/old', '--', 'python', 'review.py']
    result = command_for_host(command, 'different-node')
    assert result[0] == 'ssh' and result[-2] == 'different-node'
    args = shlex.split(result[-1])
    assert 'OMP_NUM_THREADS=1' in args and 'OPENBLAS_NUM_THREADS=1' in args
    program = args[-1]
    ast.parse(program)
    assert 'os.sched_setaffinity' in program
    assert repr(command) in program
    assert 'gpu_reserve' not in program
