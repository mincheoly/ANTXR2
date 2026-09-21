"""Run antxr2_coexpr.py over many populations as concurrent subprocesses.

One process per population (num_cpus=1 inside memento) -- on Windows that beats
memento's own multiprocessing, whose spawn overhead ate most of the speedup in
benchmarking. Concurrency is capped both by worker count and by a budget on the
total nnz in flight, because peak RSS tracks matrix size. Resumable: a
population whose .csv already exists is skipped.

Usage:
    python antxr2_coexpr_driver.py <indir> <outdir> [--workers 5]
        [--nnz-budget 120e6] [--boot 2000] [--q 0.15] [--pattern ""]
"""
import argparse
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RUNNER = os.path.join(HERE, "antxr2_coexpr.py")


def child_env():
    """Children need the conda env's DLL directories on PATH: launched without
    them the worker dies silently at a delay-loaded DLL (exit 0xC06D007E) with
    no traceback, even though a bare `import memento` succeeds."""
    root = os.path.dirname(os.path.abspath(sys.executable))
    dirs = [root, os.path.join(root, "Library", "bin"),
            os.path.join(root, "Library", "mingw-w64", "bin"),
            os.path.join(root, "Library", "usr", "bin"),
            os.path.join(root, "Scripts")]
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join(dirs + [env.get("PATH", "")])
    return env


def pop_nnz(path):
    z = np.load(path, allow_pickle=False)
    return int(z["indptr"][-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("outdir")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--nnz-budget", type=float, default=120e6)
    ap.add_argument("--boot", type=int, default=2000)
    ap.add_argument("--q", type=float, default=0.15)
    ap.add_argument("--pattern", default="")
    ap.add_argument("--max-celltypes", type=int, default=None)
    ap.add_argument("--skip", default="", help="comma-separated populations to leave alone")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    global CENV
    CENV = child_env()

    skip = {s for s in a.skip.split(",") if s}
    todo = []
    for fn in sorted(os.listdir(a.indir)):
        if not fn.endswith(".npz") or a.pattern not in fn:
            continue
        pop = fn[:-4]
        # .claim lets a second driver instance share the same output directory
        if pop in skip or any(os.path.exists(os.path.join(a.outdir, f"{pop}.{e}"))
                              for e in ("csv", "json", "claim")):
            continue
        p = os.path.join(a.indir, fn)
        nnz = pop_nnz(p)
        todo.append((nnz, pop, p))
    todo.sort()                      # small first: fast feedback, tail is the big ones
    print(f"{len(todo)} populations queued, total nnz {sum(t[0] for t in todo)/1e6:.0f}M",
          flush=True)

    running = []                     # (proc, nnz, pop, t0)
    t_start = time.time()
    done = 0
    while todo or running:
        while todo:
            nnz, pop, p = todo[0]
            live = sum(r[1] for r in running)
            if running and (len(running) >= a.workers or live + nnz > a.nnz_budget):
                break
            todo.pop(0)
            # claim atomically at dispatch time, not at queue-build time: a
            # concurrent driver instance may have taken this population since
            try:
                os.close(os.open(os.path.join(a.outdir, f"{pop}.claim"),
                                 os.O_CREAT | os.O_EXCL | os.O_WRONLY))
            except FileExistsError:
                continue
            if os.path.exists(os.path.join(a.outdir, f"{pop}.csv")):
                continue
            cmd = [sys.executable, RUNNER, p, a.outdir, "--boot", str(a.boot),
                   "--cpus", "1", "--q", str(a.q)]
            if a.max_celltypes:
                cmd += ["--max-celltypes", str(a.max_celltypes)]
            log = open(os.path.join(a.outdir, f"{pop}.log"), "w")
            proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=CENV)
            running.append((proc, nnz, pop, time.time(), log))
            print(f"[{time.time()-t_start:7.0f}s] START {pop} (nnz {nnz/1e6:.1f}M, "
                  f"{len(running)} running)", flush=True)
        time.sleep(5)
        still = []
        for proc, nnz, pop, t0, log in running:
            if proc.poll() is None:
                still.append((proc, nnz, pop, t0, log))
            else:
                log.close()
                done += 1
                print(f"[{time.time()-t_start:7.0f}s] {'DONE' if proc.returncode == 0 else 'FAIL'} "
                      f"{pop} rc={proc.returncode} {(time.time()-t0)/60:.1f} min "
                      f"({done} finished, {len(todo)} queued)", flush=True)
        running = still
    print(f"ALL DONE in {(time.time()-t_start)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
