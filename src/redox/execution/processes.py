"""Release state-worker children without leaving reusable scientific clients stale."""
import sys

import psutil


def stop_owned_children(timeout=10):
    """Close the CPU thermodynamics client, then stop remaining owned children.

    GPU queue processes reuse Python modules between molecular states. Killing
    the cached thermodynamic subprocess without closing its client leaves the
    next state writing to a dead pipe. No thermodynamic code or protocol changes
    are needed: the client's existing close routine resets that ownership.

    Call only at a state boundary or after interrupting its computation. This
    function acts on children of this process; it never releases GPU leases.
    """
    if timeout<=0:raise ValueError('positive child shutdown timeout required')
    thermal=sys.modules.get('redox.qm.qrrho')
    if thermal is not None:
        thermal._close_worker()
    children=psutil.Process().children(recursive=True)
    for child in reversed(children):
        try:child.terminate()
        except psutil.NoSuchProcess:pass
    _,alive=psutil.wait_procs(children,timeout=timeout)
    for child in alive:
        try:child.kill()
        except psutil.NoSuchProcess:pass
    _,survivors=psutil.wait_procs(alive,timeout=timeout)
    if survivors:raise RuntimeError('Owned state-worker children did not exit')
