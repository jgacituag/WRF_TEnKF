#!/usr/bin/env python3
"""Did da_update_bc read the analysis? Check a wrfutil DA cycle step by step.

At each forecast step, run_ensfcst.sh renames real.exe's wrfinput_d01 to
wrfinput_d01.org, copies the analysis to wrfinput_d01 and runs da_update_bc.
da_update_bc reads the analysis (da_file) to update the lateral boundaries and
writes the lower-boundary and land-surface fields into it. On Hydra the rename
and the copy ran on the job's head node, while da_update_bc ran on the member's
node. That node could still resolve the name wrfinput_d01 to real.exe's file
(November 2026). da_update_bc then built the boundary conditions from real.exe's
fields and wrote the surface update into wrfinput_d01.org, partly or fully.
wrf.exe still started from the analysis.

da_update_bc's log prints the first element of every field it reads
("Sample data="). For each step s and member, this script compares the samples
of the state fields read from ./wrfinput_d01 (HIST/GUES/<t_s>/
da_update_bc_<s>_<mem>.log) with the first element of the same fields in the
analysis the step started from (HIST/ANAL/<t_s-1>/anal<mem>), and with the
guess valid then (HIST/GUES/<t_s-1>/gues<mem>):
  ok       every sampled state field equals the analysis
  WRONG    some sampled field differs from the analysis: da_update_bc read
           another file. Reported as k/n fields; "g" counts the wrong samples
           that equal the guess, which carries real.exe's boundary values at
           the domain corner.
  missing  log or analysis not found
A wrong read is invisible only where real.exe's file and the analysis agree at
the domain corner, e.g. an analysis that leaves the corner unchanged and
boundaries made from that same analysis.

Usage:
  check_anal_reads.py EXP [--steps S0-S1] [--src DIR ...] [--jobs N] [--members]
EXP is an experiment folder, on scratch or archived, with HIST/ (and conf/exp.conf
for ANALYSIS_FREQ, default 300 s). --src adds experiment folders searched for
HIST/ANAL/<t> and HIST/GUES/<t> that EXP lacks, e.g. the source of a restart (its
first analysis, and the guess at that time). --members lists the members per step.
Exit status: 0 if no WRONG member, 1 if any, 2 if nothing was checked.
Runs on the login node: Python 3 with netCDF4 (~/salidas/miniconda3/bin/python3).
The fix: bin/run_ensfcst.sh (WRF_TEnKF) renames and copies on the member's node.
"""
import argparse
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timedelta
from pathlib import Path

import netCDF4

STATE = ("MU", "U", "V", "T", "PH", "QVAPOR")   # fields da_update_bc reads from da_file
RTOL = 1e-9   # the log prints 12 significant digits of a float32: exact to ~1e-12


def log_samples(path):
    """First 'Sample data' of each state field read from ./wrfinput_d01."""
    out, cur_file, cur_var = {}, None, None
    with open(path, errors="replace") as f:
        for line in f:
            if "open netcdf file" in line:
                cur_file = line.split()[-1]
            elif "varnam=" in line:
                cur_var = line.split("varnam=")[1].strip()
            elif "Sample data=" in line and cur_file == "./wrfinput_d01" \
                    and cur_var in STATE and cur_var not in out:
                out[cur_var] = float(line.split("=")[1])
    return out


def first_elements(path, names):
    with netCDF4.Dataset(path) as d:
        d.set_auto_mask(False)
        return {v: float(d[v][(0,) * d[v].ndim]) for v in names if v in d.variables}


def same(x, ref):
    return abs(x - ref) <= RTOL * max(abs(ref), 1e-30)


def find(roots, kind, tag, name):
    for r in roots:
        p = r / "HIST" / kind / tag / name
        if p.exists():
            return p
    return None


def check_member(args):
    step, mem, log, anal, gues = args
    if log is None or anal is None:
        return step, mem, "missing", 0, 0, 0
    s = log_samples(log)
    if not s:
        return step, mem, "missing", 0, 0, 0
    a = first_elements(anal, s)
    g = first_elements(gues, s) if gues is not None else {}
    wrong = [v for v in s if v in a and not same(s[v], a[v])]
    as_guess = [v for v in wrong if v in g and same(s[v], g[v])]
    return step, mem, ("WRONG" if wrong else "ok"), len(wrong), len(s), len(as_guess)


def conf_value(exp, name, default):
    f = exp / "conf" / "exp.conf"
    if f.exists():
        for line in f.read_text(errors="replace").splitlines():
            m = re.match(rf"\s*(export\s+)?{name}=['\"]?([^'\"#\s]+)", line)
            if m:
                return m.group(2)
    return default


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("exp", type=Path)
    ap.add_argument("--steps", help="S0-S1 (default: every step with logs)")
    ap.add_argument("--src", type=Path, action="append", default=[])
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--members", action="store_true")
    o = ap.parse_args()
    exp = o.exp.resolve()
    roots = [exp] + [p.resolve() for p in o.src]
    freq = int(conf_value(exp, "ANALYSIS_FREQ", "300"))
    s0, s1 = (int(x) for x in o.steps.split("-")) if o.steps else (0, 10 ** 6)

    tasks = []
    for tdir in sorted((exp / "HIST" / "GUES").glob("[0-9]" * 14)):
        for log in sorted(tdir.glob("da_update_bc_*_*.log")):
            m = re.match(r"da_update_bc_(\d+)_(\d+)\.log$", log.name)
            if not m or not s0 <= int(m.group(1)) <= s1:
                continue
            step, mem = int(m.group(1)), int(m.group(2))
            t1 = datetime.strptime(tdir.name, "%Y%m%d%H%M%S")
            t0 = (t1 - timedelta(seconds=freq)).strftime("%Y%m%d%H%M%S")
            tasks.append((step, mem, log, find(roots, "ANAL", t0, f"anal{mem:05d}"),
                          find(roots, "GUES", t0, f"gues{mem:05d}")))
    if not tasks:
        print(f"No da_update_bc logs in {exp}/HIST/GUES for steps {o.steps or 'all'}")
        return 2

    with ProcessPoolExecutor(max(1, o.jobs)) as pool:
        res = sorted(pool.map(check_member, tasks))

    print(f"{exp}  (state fields {', '.join(STATE)}; WRONG = da_update_bc did not read the analysis)")
    print(f"{'step':>4} {'from':>5} {'members':>7} {'ok':>4} {'WRONG':>5} {'missing':>7}")
    tot = {"ok": 0, "WRONG": 0, "missing": 0}
    for step in sorted({r[0] for r in res}):
        rs = [r for r in res if r[0] == step]
        c = {k: sum(r[2] == k for r in rs) for k in tot}
        for k in tot:
            tot[k] += c[k]
        t1 = next(t[2].parent.name for t in tasks if t[0] == step)
        t0 = datetime.strptime(t1, "%Y%m%d%H%M%S") - timedelta(seconds=freq)
        line = f"{step:>4} {t0:%H:%M} {len(rs):>7} {c['ok']:>4} {c['WRONG']:>5} {c['missing']:>7}"
        if o.members and c["WRONG"]:
            line += "   " + " ".join(f"{r[1]:03d}({r[3]}/{r[4]},g{r[5]})" for r in rs if r[2] == "WRONG")
        print(line)
    n = sum(tot.values())
    print(f"total: {n} member-steps, {tot['ok']} ok, {tot['WRONG']} WRONG, {tot['missing']} missing")
    return 1 if tot["WRONG"] else 0


if __name__ == "__main__":
    sys.exit(main())
