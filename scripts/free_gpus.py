#!/usr/bin/env python
"""Pick GPUs that are ACTUALLY idle — across this SHARED node and (optionally) peer nodes.

The lambda cluster (lambda1..lambda4, and this node lambda2) shares one NFS filesystem and
one conda env, so a job launched on any reachable node reads/writes the SAME paths. But the
node is SHARED with other users: never introduce contention (see CLAUDE.md). A GPU counts as
free only if utilization <= UTIL_MAX (%) AND used memory <= MEM_MAX (MiB).

  python scripts/free_gpus.py                         # local idle indices -> "2,5"
  python scripts/free_gpus.py -n 1                     # at most 1        -> "2"
  python scripts/free_gpus.py --list                  # local human table
  python scripts/free_gpus.py --hosts lambda2,lambda4 # cross-node -> lines "host idx"
  python scripts/free_gpus.py --hosts lambda1,lambda2,lambda4 --list

Cross-node output is one "host idx" per line (unreachable hosts are skipped with a note on
stderr), so a launcher can round-robin jobs onto free (host, gpu) slots.
"""
from __future__ import annotations
import argparse
import subprocess
import sys

UTIL_MAX = 5      # percent
MEM_MAX = 500     # MiB

_QUERY = ["nvidia-smi", "--query-gpu=index,utilization.gpu,memory.used",
          "--format=csv,noheader,nounits"]


def _parse(text):
    rows = []
    for line in text.strip().splitlines():
        parts = [x.strip() for x in line.split(",")]
        if len(parts) != 3 or not parts[0].isdigit():
            continue
        rows.append((int(parts[0]), int(parts[1]), int(parts[2])))
    return rows


def query_local():
    try:
        return _parse(subprocess.check_output(_QUERY, text=True))
    except (FileNotFoundError, subprocess.CalledProcessError):
        return []


def query_remote(host, timeout=15):
    try:
        out = subprocess.check_output(
            ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=6",
             "-o", "StrictHostKeyChecking=no", host, " ".join(_QUERY)],
            text=True, stderr=subprocess.DEVNULL, timeout=timeout)
        return _parse(out)
    except Exception:
        return None   # unreachable / no GPUs


def free(rows):
    return [idx for idx, util, mem in rows if util <= UTIL_MAX and mem <= MEM_MAX]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=None, help="cap the number of GPUs returned")
    ap.add_argument("--list", action="store_true", help="human-readable table")
    ap.add_argument("--hosts", default=None,
                    help="comma-sep hosts to scan over SSH (e.g. lambda1,lambda2,lambda4); "
                         "omit for local-only")
    args = ap.parse_args()

    import socket
    me = socket.gethostname().split(".")[0]

    if not args.hosts:
        rows = query_local()
        if args.list:
            for idx, util, mem in rows:
                tag = "FREE" if (util <= UTIL_MAX and mem <= MEM_MAX) else "busy"
                print(f"  GPU{idx}: util={util:3d}%  mem={mem:6d} MiB  [{tag}]")
            return
        idle = free(rows)
        if args.n is not None:
            idle = idle[:args.n]
        print(",".join(str(i) for i in idle))
        sys.exit(0 if idle else 1)

    # cross-node
    slots = []   # (host, idx)
    for host in [h.strip() for h in args.hosts.split(",") if h.strip()]:
        rows = query_local() if host == me else query_remote(host)
        if rows is None:
            print(f"[skip] {host}: unreachable", file=sys.stderr)
            continue
        hfree = free(rows)
        if args.list:
            print(f"{host}: free {hfree or '(none)'}")
        for idx in hfree:
            slots.append((host, idx))
    if args.list:
        return
    if args.n is not None:
        slots = slots[:args.n]
    for host, idx in slots:
        print(f"{host} {idx}")
    sys.exit(0 if slots else 1)


if __name__ == "__main__":
    main()
