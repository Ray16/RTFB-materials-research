#!/usr/bin/env python
"""Pick GPUs that are ACTUALLY free — util/mem low AND no OTHER user resident AND host allowed.

The lambda cluster shares one NFS filesystem and one conda env, so a job launched on any reachable
node reads/writes the SAME paths. But nodes are SHARED with other users, so a GPU counts as free
ONLY if ALL of:
  * utilization <= UTIL_MAX (%)              (cluster.env)
  * used memory <= MEM_MAX (MiB)             (cluster.env)
  * NO process owned by another user is resident on it   (ownership check — the key fix)
  * its host is in ALLOWED_HOSTS and NOT in RESERVED_HOSTS   (cluster.env; e.g. lambda13 reserved)

Ownership is decided by scripts/fleet/gpu_probe.sh (idle-but-resident foreign jobs count as OCCUPIED
— an idle 314 MiB foreign context is still someone else's GPU). Policy lives in scripts/fleet/cluster.env.

  python scripts/free_gpus.py                         # local free indices -> "2,5"
  python scripts/free_gpus.py -n 1                     # at most 1        -> "2"
  python scripts/free_gpus.py --list                  # local human table
  python scripts/free_gpus.py --hosts lambda1,lambda2 # cross-node -> lines "host idx"
  python scripts/free_gpus.py --all                    # scan every ALLOWED host in cluster.env
  python scripts/free_gpus.py --all --list
"""
from __future__ import annotations
import argparse
import os
import socket
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROBE = os.path.join(ROOT, "scripts", "fleet", "gpu_probe.sh")
CLUSTER_ENV = os.path.join(ROOT, "scripts", "fleet", "cluster.env")


def load_env():
    """Parse scripts/fleet/cluster.env (plain KEY="value") with sane fallbacks."""
    cfg = {"ALLOWED_HOSTS": "", "RESERVED_HOSTS": "", "UTIL_MAX": "5", "MEM_MAX": "500"}
    try:
        with open(CLUSTER_ENV) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#") or "=" not in line:
                    continue
                k, v = line.split("=", 1)
                cfg[k.strip()] = v.strip().strip('"').strip("'")
    except FileNotFoundError:
        pass
    return cfg


CFG = load_env()
UTIL_MAX = int(CFG["UTIL_MAX"])
MEM_MAX = int(CFG["MEM_MAX"])
ALLOWED = CFG["ALLOWED_HOSTS"].split()
RESERVED = set(CFG["RESERVED_HOSTS"].split())


def _parse(text):
    """gpu_probe.sh lines: idx,util,mem,foreign_owners -> (idx, util, mem, foreign_str)."""
    rows = []
    for line in text.strip().splitlines():
        parts = line.split(",", 3)
        if len(parts) < 4 or not parts[0].strip().isdigit():
            continue
        rows.append((int(parts[0]), int(parts[1]), int(parts[2]), parts[3].strip()))
    return rows


def probe_local():
    try:
        return _parse(subprocess.check_output(["bash", PROBE], text=True))
    except (FileNotFoundError, subprocess.CalledProcessError):
        return []


def probe_remote(host, timeout=20):
    try:
        with open(PROBE) as fh:
            out = subprocess.check_output(
                ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=6",
                 "-o", "StrictHostKeyChecking=no", host, "bash -s"],
                stdin=fh, text=True, stderr=subprocess.DEVNULL, timeout=timeout)
        return _parse(out)
    except Exception:
        return None   # unreachable / no GPUs


def free(rows):
    """A GPU is free only if quiet AND owned by nobody else."""
    return [idx for idx, util, mem, foreign in rows
            if util <= UTIL_MAX and mem <= MEM_MAX and not foreign]


def host_allowed(host):
    short = host.split(".")[0]
    if short in RESERVED:
        return False
    return (not ALLOWED) or (short in ALLOWED)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=None, help="cap the number of GPUs/slots returned")
    ap.add_argument("--list", action="store_true", help="human-readable table")
    ap.add_argument("--hosts", default=None,
                    help="comma-sep hosts to scan over SSH; omit for local-only")
    ap.add_argument("--all", action="store_true",
                    help="scan every ALLOWED host in cluster.env")
    args = ap.parse_args()

    me = socket.gethostname().split(".")[0]

    # local-only mode
    if not args.hosts and not args.all:
        rows = probe_local()
        if args.list:
            for idx, util, mem, foreign in rows:
                tag = "FREE" if (util <= UTIL_MAX and mem <= MEM_MAX and not foreign) else "busy"
                extra = f"  other-user={foreign}" if foreign else ""
                print(f"  GPU{idx}: util={util:3d}%  mem={mem:6d} MiB  [{tag}]{extra}")
            return
        idle = free(rows)
        if args.n is not None:
            idle = idle[:args.n]
        print(",".join(str(i) for i in idle))
        sys.exit(0 if idle else 1)

    # cross-node mode
    hosts = ALLOWED if args.all else [h.strip() for h in args.hosts.split(",") if h.strip()]
    slots = []   # (host, idx)
    for host in hosts:
        short = host.split(".")[0]
        if short in RESERVED:
            print(f"[skip] {host}: RESERVED (never launch here)", file=sys.stderr)
            continue
        if ALLOWED and short not in ALLOWED:
            print(f"[skip] {host}: not in ALLOWED_HOSTS", file=sys.stderr)
            continue
        rows = probe_local() if short == me else probe_remote(host)
        if rows is None:
            print(f"[skip] {host}: unreachable", file=sys.stderr)
            continue
        hfree = free(rows)
        if args.list:
            occ = {idx: foreign for idx, _u, _m, foreign in rows if foreign}
            note = f"  (other users: {occ})" if occ else ""
            print(f"{host}: free {hfree or '(none)'}{note}")
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
