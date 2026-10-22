#!/usr/bin/env python
"""Source-preserving launcher fix for an already running frozen native campaign.

The frozen module's compound `cd && nohup ... &` leaves a waiting subshell
attached to SSH. Its worker uses absolute paths and the shared GPU wrapper
sets the compute cwd, so launch only the redirected nohup command. Numerical
producer bytes remain unchanged. Retire this adapter after the campaign closes
and incorporate the same fix in the reusable launch function.
"""

import argparse
from contextlib import contextmanager
import shlex

from redox.core.common import ROOT
from redox.workflows import benchmark_native


def detached_argv(argv):
    prefix = "cd " + shlex.quote(str(ROOT)) + " && "
    if (
        argv
        and argv[0] == "ssh"
        and argv[-1].startswith(prefix + "nohup ")
        and "redox.workflows.benchmark_native worker" in argv[-1]
    ):
        return [*argv[:-1], argv[-1][len(prefix) :]]
    return argv


@contextmanager
def corrected_launch():
    original = benchmark_native.subprocess.run

    def run(argv, **kwargs):
        return original(detached_argv(argv), **kwargs)

    benchmark_native.subprocess.run = run
    try:
        yield
    finally:
        benchmark_native.subprocess.run = original


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True)
    args = parser.parse_args()
    with corrected_launch():
        benchmark_native.launch(args.campaign)
