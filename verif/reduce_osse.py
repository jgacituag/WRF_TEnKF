"""Reduce one OSSE run (one ensemble) to the sums and compact fields the verification notebooks need.

One call reduces ONE run at the requested valid times, reading only that run's members and
the truth, so that each forecast run can be reduced as soon as it ends and then leave the
scratch area. Everything that compares several runs (the echo domain shared by all of them,
false and missed echo, the fraction of cells that improve, the maps) is computed in the
notebooks from the compact fields (vcommon.py).

Runs and valid times (DATE = 2024-03-19, analysis at 18:00):
  --exp noda   18:00        "prior": DA run guesses  <osse>/da4km/HIST/GUES/<DATE>180000/gues<NNNNN>
               18:05-19:00  forecast without DA      <osse>/osse18/fc_noda/HIST/DAFCST/<DATE>180000/<MMM>/wrfout_d01_<time>
  --exp ntN    18:00        OSSE analyses            <osse>/osse18/anal_ntN/anal<NNNNN>
               18:05-19:00  forecast from them       <osse>/osse18/fc_ntN/... (same layout)
  truth        member TRUTH: its guess at 18:00, then member TRUTH of fc_noda.
At 18:00 the run "noda" is written under the name "prior".

Output files, in --out (default <osse>/osse18/verif):
  stats_<exp>_<HHMM>.csv    pooled sums per (variable, domain): n, sum of squared error of the
                            ensemble mean (sse), sum of error (se), sum of ensemble variance
                            with ddof = 1 (svar) and sum of CRPS (scrps). Ratios are formed
                            after pooling, in the notebook (vcommon.load_stats).
  levels_<exp>_<HHMM>.npz   the same sums per model level, key "<var>|<domain>", array
                            (5, nk) in the order n, sse, se, svar, scrps.
  fields_<exp>_<HHMM>.npz   per-cell fields for the multi-run diagnostics, on the box:
                            ref_mean, ref_var, ref_crps   reflectivity [dBZ], (nk, ny, nx) float32
                            hyd_mean                      ensemble mean of qr + qs + qg [kg/kg]
                            colmax                        column-maximum reflectivity of every
                                                          member, (NE, ny, nx) float16
                            colmax_mean                   column maximum of ref_mean, (ny, nx)
  truth_<HHMM>.npz          the truth on the box: ref, hyd (qr + qs + qg), colmax, storm
                            (columns with colmax >= STORM dBZ). Written once per time.
  masks.npz                 the verification box (box, region), the observed cells, mean
                            level heights (zmean), lev5, lat, lon, n_obs and jsec (the row
                            with the most storm columns at 18:00, for vertical sections).

Domains of the sums (each fixed by the truth or the box alone, so the same for every run):
  global      every cell of the verification box: interior points (3..n-2 in Fortran
              numbering), levels 2..nlev-2 with height <= ZMAX, as the synthetic observations
  observed    box cells carrying a synthetic observation (rad0101.dat of the run --obs-from)
  storm       box cells in columns where the truth's column maximum reaches STORM dBZ
  truth_echo  box cells where the truth reaches ECHO dBZ

Reflectivity is computed member by member with the Tong & Xue operator of the LETKF
(osse_io.calc_ref; check it with check_operator.py) and averaged in dBZ (mean of H(x_i)).

Usage (one node: on Hydra osse/jobs/reduce_osse_hydra.sh):
  python3 verif/reduce_osse.py --osse ~/salidas/osse --exp nt1 --times 1800
  python3 verif/reduce_osse.py --osse ~/salidas/osse --exp noda          # 18:00 to 19:00
  options: --out DIR  --ne 59 --truth 60 --zmax 16000 --echo 10 --storm 20 --obs-from nt1
           --nproc 16 --mem-frac 0.5 --kchunk 0 --timeout 900 --overwrite
Memory: the levels are processed in chunks sized to the memory limit (cgroup, else the
node's memory), so one chunk of the ensemble takes at most --mem-frac of it. If a reading
process dies (memory), the script stops after --timeout s instead of waiting forever.
"""
import argparse
import csv
import os
import sys
import time
from datetime import datetime
from multiprocessing import TimeoutError as PoolTimeout, get_context
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import osse_io as io  # noqa: E402

DOMAINS = ("global", "observed", "storm", "truth_echo")
SUM_KEYS = ("n", "sse", "se", "svar", "scrps")
HYD = ("qr", "qs", "qg")


# ----------------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------------
class Layout:
    def __init__(self, osse, date, ne, truth):
        self.osse, self.date, self.ne, self.truth = Path(osse), date, ne, truth
        self.t0 = datetime.strptime(date + "1800", "%Y%m%d%H%M")

    def gues(self, m):
        return self.osse / "da4km/HIST/GUES" / f"{self.date}180000" / f"gues{m:05d}"

    def anal(self, exp, m):
        return self.osse / "osse18" / f"anal_{exp}" / f"anal{m:05d}"

    def fcst(self, exp, m, t):
        return (self.osse / "osse18" / f"fc_{exp}" / "HIST/DAFCST" / f"{self.date}180000"
                / f"{m:03d}" / f"wrfout_d01_{t:%Y-%m-%d_%H:%M:%S}")

    def files(self, exp, t):
        """Member files of run `exp` ('noda' or an OSSE run 'ntN') valid at t."""
        mem = range(1, self.ne + 1)
        if t == self.t0:
            return [self.gues(m) if exp == "noda" else self.anal(exp, m) for m in mem]
        return [self.fcst(exp, m, t) for m in mem]

    def truth_file(self, t):
        return self.gues(self.truth) if t == self.t0 else self.fcst("noda", self.truth, t)

    def name(self, exp, t):
        """The run's name in the output files: the reference at 18:00 is the prior."""
        return "prior" if exp == "noda" and t == self.t0 else exp


# ----------------------------------------------------------------------------------
# Reading an ensemble, one worker process per member file
# ----------------------------------------------------------------------------------
def _read_member(args):
    m, path, box = args
    return m, io.read_wrf(path, box)


def mem_limit_gb():
    """Memory this process may use: the cgroup limit of the session, or the physical
    memory when there is none."""
    for f in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            v = open(f).read().strip()
            if v != "max" and int(v) < 1 << 60:
                return int(v) / 2 ** 30
        except (OSError, ValueError):
            pass
    return os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 2 ** 30


class Reader:
    """A pool of worker processes (spawned, not forked: HDF5 state is not inherited)
    that read member files and return their arrays."""

    def __init__(self, nproc, timeout):
        self.pool = get_context("spawn").Pool(nproc)
        self.timeout = timeout

    def ensemble(self, files, box, label=""):
        """(NE, nk, ny, nx) float32 arrays per variable."""
        missing = [str(f) for f in files if not Path(f).exists()]
        if missing:
            raise FileNotFoundError(f"{len(missing)} member files missing, first {missing[0]}")
        k0, k1, j0, j1, i0, i1 = box
        shape = (len(files), k1 - k0, j1 - j0, i1 - i0)
        x = {v: np.empty(shape, np.float32) for v in io.VARS}
        it = self.pool.imap_unordered(_read_member, [(m, f, box) for m, f in enumerate(files)])
        for _ in range(len(files)):
            try:
                m, xm = it.next(timeout=self.timeout)
            except PoolTimeout as err:
                raise RuntimeError(f"reading {label}: no member arrived in {self.timeout} s; "
                                   f"a reading process may have been killed for memory: "
                                   f"lower --mem-frac or --nproc") from err
            for v in io.VARS:
                x[v][m] = xm[v]
            del xm
        return x

    def close(self):
        self.pool.terminate()


# ----------------------------------------------------------------------------------
# Ensemble statistics (Chapter 4 definitions)
# ----------------------------------------------------------------------------------
def crps_sorted(ens, truth):
    """Ensemble CRPS, mean|x_i - y| - 0.5 mean|x_i - x_j|, via sorted members
    (da.metrics.crps_ensemble_sorted), member axis 0, float64 accumulation."""
    ne = ens.shape[0]
    term1 = np.zeros(truth.shape, np.float64)
    for i in range(ne):
        term1 += np.abs(ens[i] - truth)
    term1 /= ne
    s = np.sort(ens, axis=0)
    term2 = np.zeros(truth.shape, np.float64)
    for i in range(ne):
        term2 += (2 * (i + 1) - ne - 1) * s[i].astype(np.float64)
    del s
    return (term1 - term2 / ne ** 2).astype(np.float32)


def field_stats(ens, truth):
    mean = ens.mean(axis=0, dtype=np.float64)
    var = ens.var(axis=0, ddof=1, dtype=np.float64).astype(np.float32)
    err = (mean - truth).astype(np.float32)
    return dict(mean=mean.astype(np.float32), var=var, err=err, crps=crps_sorted(ens, truth))


def level_sums(st, dom):
    """Pooled sums of one variable over one domain, per level: dict key -> (nk,)."""
    w = dom.astype(np.float64)
    e = st["err"].astype(np.float64)
    return dict(n=w.sum(axis=(1, 2)), sse=(w * e * e).sum(axis=(1, 2)), se=(w * e).sum(axis=(1, 2)),
                svar=(w * st["var"]).sum(axis=(1, 2)), scrps=(w * st["crps"]).sum(axis=(1, 2)))


def save_npz_atomic(path, **arrays):
    """Write path through a temporary file, so that a concurrent reader never sees half a file."""
    tmp = path.with_name(path.name + f".tmp{os.getpid()}.npz")
    np.savez_compressed(tmp, **arrays)
    os.replace(tmp, path)


# ----------------------------------------------------------------------------------
# Masks and truth
# ----------------------------------------------------------------------------------
def make_masks(L, a, log):
    info = io.grid_info(L.truth_file(L.t0))
    nlev, nlat, nlon = info["nlev"], info["nlat"], info["nlon"]
    zw = info["zw"][:nlev]
    kmax = max(k for k in range(1, nlev - 2) if zw[k, 2:nlat - 2, 2:nlon - 2].min() <= a.zmax)
    box = (1, kmax + 1, 2, nlat - 2, 2, nlon - 2)
    k0, k1, j0, j1, i0, i1 = box
    region = zw[k0:k1, j0:j1, i0:i1] <= a.zmax
    obs = io.read_radar_obs(L.osse / "osse18" / f"anal_{a.obs_from}" / "rad0101.dat")
    kk, jj, ii, _ = io.obs_to_ijk(obs, info)
    inbox = (kk >= k0) & (kk < k1) & (jj >= j0) & (jj < j1) & (ii >= i0) & (ii < i1)
    if not inbox.all():
        log(f"{(~inbox).sum()} observations outside the verification box, not in 'observed'")
    observed = np.zeros(region.shape, bool)
    observed[kk[inbox] - k0, jj[inbox] - j0, ii[inbox] - i0] = True
    zmean = zw[k0:k1, j0:j1, i0:i1].mean(axis=(1, 2))
    lev5 = int(np.argmin(np.abs(zmean - 5000.0)))
    truth = io.read_wrf(L.truth_file(L.t0), box)
    tcol = np.where(region, truth["ref"], -np.inf).max(axis=0)
    storm = np.broadcast_to(tcol >= a.storm, region.shape) & region
    jsec = int(np.argmax(storm[lev5].sum(axis=1)))
    log(f"box k {k0}:{k1} j {j0}:{j1} i {i0}:{i1}; {region.sum()} cells, "
        f"{observed.sum()} observed ({len(obs)} observations); section row {jsec}")
    return dict(box=np.array(box), region=region, observed=observed & region, zmean=zmean,
                lev5=lev5, jsec=jsec, lat=info["lat"][j0:j1, i0:i1], lon=info["lon"][j0:j1, i0:i1],
                n_obs=len(obs), zmax=a.zmax, echo=a.echo, storm=a.storm)


def truth_fields(truth, region, storm_dbz):
    ref = np.where(region, truth["ref"], 0.0).astype(np.float32)
    colmax = ref.max(axis=0)
    hyd = np.where(region, truth["qr"] + truth["qs"] + truth["qg"], 0.0).astype(np.float32)
    return dict(ref=ref, hyd=hyd, colmax=colmax, storm=colmax >= storm_dbz)


# ----------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--osse", required=True, help="folder holding da4km/ and osse18/")
    ap.add_argument("--exp", required=True, help="run: noda or ntN")
    ap.add_argument("--date", default="20240319")
    ap.add_argument("--times", default=",".join(f"{18 + m // 60:02d}{m % 60:02d}" for m in range(0, 61, 5)))
    ap.add_argument("--out", default=None)
    ap.add_argument("--ne", type=int, default=59)
    ap.add_argument("--truth", type=int, default=60)
    ap.add_argument("--zmax", type=float, default=16000.0)
    ap.add_argument("--echo", type=float, default=10.0, help="threshold of the truth_echo domain [dBZ]")
    ap.add_argument("--storm", type=float, default=20.0, help="truth column maximum of the storm domain [dBZ]")
    ap.add_argument("--obs-from", default="nt1", help="OSSE run whose rad0101.dat gives the observed cells")
    ap.add_argument("--nproc", type=int, default=16)
    ap.add_argument("--mem-frac", type=float, default=0.5,
                    help="fraction of the memory one ensemble chunk may take")
    ap.add_argument("--kchunk", type=int, default=0, help="levels per chunk (0: from the memory)")
    ap.add_argument("--timeout", type=float, default=900.0, help="seconds to wait for one member")
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    if a.exp != "noda" and not (a.exp.startswith("nt") and a.exp[2:].isdigit()):
        sys.exit(f"--exp {a.exp}: expected noda or ntN")

    L = Layout(a.osse, a.date, a.ne, a.truth)
    out = Path(a.out) if a.out else L.osse / "osse18" / "verif"
    out.mkdir(parents=True, exist_ok=True)
    times = [datetime.strptime(a.date + hm, "%Y%m%d%H%M") for hm in a.times.split(",")]
    tick = time.time()
    log = lambda msg: print(f"[{time.time() - tick:7.1f} s] {a.exp} {msg}", flush=True)  # noqa: E731

    # ---- verification box, observed cells, section row (once for all runs) ----------
    mfile = out / "masks.npz"
    if mfile.exists() and not a.overwrite:
        mk = dict(np.load(mfile))
        log(f"masks from {mfile}")
    else:
        mk = make_masks(L, a, log)
        save_npz_atomic(mfile, **mk)
    box = tuple(int(b) for b in mk["box"])
    k0, k1, j0, j1, i0, i1 = box
    nk = k1 - k0
    region = mk["region"]

    # ---- level chunks that fit in memory -------------------------------------------
    gb_level = len(io.VARS) * a.ne * (j1 - j0) * (i1 - i0) * 4 / 2 ** 30
    mem = mem_limit_gb()
    cores = len(os.sched_getaffinity(0))
    if a.nproc > cores:
        log(f"--nproc {a.nproc} > {cores} cores: using {cores}")
        a.nproc = cores
    kc = a.kchunk or max(1, min(nk, int(a.mem_frac * mem / (1.6 * gb_level))))
    chunks = [(k, min(k + kc, nk)) for k in range(0, nk, kc)]
    log(f"memory {mem:.0f} GB; one level of the ensemble {gb_level:.2f} GB; "
        f"{len(chunks)} chunks of up to {kc} levels; {a.nproc} reading processes")
    reader = Reader(a.nproc, a.timeout)

    try:
        for t in times:
            tag = f"{t:%H%M}"
            name = L.name(a.exp, t)
            ffile = out / f"fields_{name}_{tag}.npz"
            if ffile.exists() and not a.overwrite:
                log(f"{tag}: {ffile.name} exists, skipped")
                continue
            files = L.files(a.exp, t)
            missing = [f for f in files if not f.exists()]
            if missing:
                raise FileNotFoundError(f"{tag}: {len(missing)} member files missing, first {missing[0]}")
            lead = int((t - L.t0).total_seconds() // 60)

            # truth of this time: read for the errors, written once
            truth = io.read_wrf(L.truth_file(t), box)
            tf = truth_fields(truth, region, a.storm)
            tfile = out / f"truth_{tag}.npz"
            if not tfile.exists() or a.overwrite:
                save_npz_atomic(tfile, **tf)
            storm3 = np.broadcast_to(tf["storm"], region.shape) & region
            domains = dict(global_=region, observed=mk["observed"], storm=storm3,
                           truth_echo=(truth["ref"] >= a.echo) & region)
            log(f"{tag} ({name}): truth read, {tf['storm'].sum()} storm columns, "
                f"{domains['truth_echo'].sum()} truth echo cells")

            acc = {}                                   # (var, domain) -> key -> (nk,)
            ref_mean = np.zeros(region.shape, np.float32)
            ref_var = np.zeros(region.shape, np.float32)
            ref_crps = np.zeros(region.shape, np.float32)
            hyd_mean = np.zeros(region.shape, np.float32)
            colmax = None
            for c in chunks:
                ens = reader.ensemble(files, (k0 + c[0], k0 + c[1], j0, j1, i0, i1), f"{name} {tag} {c}")
                sl = slice(c[0], c[1])
                for v in io.VARS:
                    st = field_stats(ens[v], truth[v][sl])
                    for dname, dom in domains.items():
                        s = level_sums(st, dom[sl])
                        d = dname.rstrip("_")
                        for key, arr in s.items():
                            acc.setdefault((v, d), {}).setdefault(key, np.zeros(nk))[sl] = arr
                    if v == "ref":
                        ref_mean[sl], ref_var[sl], ref_crps[sl] = st["mean"], st["var"], st["crps"]
                        cm = np.where(region[sl], ens["ref"], 0.0).max(axis=1)
                        colmax = cm if colmax is None else np.maximum(colmax, cm)
                    elif v in HYD:
                        hyd_mean[sl] += st["mean"]
                    del st
                del ens
                log(f"{tag}: levels {c[0]}-{c[1] - 1} done")

            for arr in (ref_mean, ref_var, ref_crps, hyd_mean):
                arr[~region] = 0.0
            rows, levels = [], {}
            for (v, d), s in acc.items():
                levels[f"{v}|{d}"] = np.stack([s[k] for k in SUM_KEYS])
                rows.append(dict(time=tag, lead=lead, exp=name, var=v, domain=d,
                                 **{k: float(np.sum(s[k])) for k in SUM_KEYS}))
            sfile = out / f"stats_{name}_{tag}.csv"
            with open(sfile.with_suffix(".tmp"), "w", newline="") as f:
                wr = csv.DictWriter(f, fieldnames=["time", "lead", "exp", "var", "domain", *SUM_KEYS])
                wr.writeheader()
                wr.writerows(rows)
            os.replace(sfile.with_suffix(".tmp"), sfile)
            save_npz_atomic(out / f"levels_{name}_{tag}.npz", **levels)
            # the fields file last: its presence marks the time as done
            save_npz_atomic(ffile, ref_mean=ref_mean, ref_var=ref_var, ref_crps=ref_crps,
                            hyd_mean=hyd_mean, colmax=colmax.astype(np.float16),
                            colmax_mean=ref_mean.max(axis=0))
            log(f"{tag}: written")
    finally:
        reader.close()


if __name__ == "__main__":
    main()
