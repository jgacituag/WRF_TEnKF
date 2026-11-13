"""Reduce the OSSE output to the numbers and slices the verification notebook needs.

For each valid time it reads the truth (member TRUTH), the reference ensemble and each
experiment's ensemble (members 1..NE), and writes:

  stats_<HHMM>.csv    pooled sums per (experiment, variable, domain): n, sum of squared
                      error of the ensemble mean, sum of error (bias), sum of ensemble
                      variance (ddof=1), sum of CRPS; and, against the reference, the
                      number of cells whose |error| improves or degrades beyond TOL.
                      The notebook forms RMSE, CRPS, spread and skill scores from them,
                      so that every ratio is taken after pooling, as in Chapter 4.
  levels_<HHMM>.npz   the same sums per model level.
  slices_<HHMM>.npz   2-D fields for figures: reflectivity of the truth and of each
                      ensemble mean at the level nearest 5 km, column maxima, a vertical
                      section, the spread at 5 km and the cell-wise skill at 5 km.
  masks.npz           the verification box and the fixed domains (from 18:00).

Times and sources (DATE = 2024-03-19, analysis at 18:00):
  18:00        reference "prior": DA run guesses   <osse>/da4km/HIST/GUES/<DATE>180000/gues<NNNNN>
               experiment ntN:    OSSE analyses    <osse>/osse18/anal_ntN/anal<NNNNN>
               truth:             guess of TRUTH in the same GUES folder
  18:05-19:00  reference "noda":  forecast without DA  <osse>/osse18/fc_noda/HIST/DAFCST/<DATE>180000/<MMM>/wrfout_d01_<time>
               experiment ntN:    <osse>/osse18/fc_ntN/... (same layout)
               truth:             member TRUTH of fc_noda

Domains (Chapter 4, nbcommon / da.metrics.domain_masks):
  global    every cell of the verification box: interior points (3..n-2 in Fortran
            numbering), levels 2..nlev-2 with height <= ZMAX, as the synthetic
            observations;
  reached   cells where the first experiment's analysis differs from the prior in any
            member and variable (the update reached them); fixed after 18:00;
  storm     columns whose truth column maximum reflectivity >= 20 dBZ, at the valid time;
  observed  cells carrying a synthetic observation (rad0101.dat of the first experiment);
            fixed after 18:00.

Reflectivity is computed member by member with the Tong & Xue operator of the LETKF
(osse_io.calc_ref; check it with check_operator.py) and averaged in dBZ.

Usage (on a Fugaku pre/post node, python3.12 with netCDF4 in PYLIB):
  python3.12 verif/reduce_osse.py --osse /data/ra000007/$USER/osse --exps nt1,nt2,nt4 --nproc 32
  options: --times 1800,1805,...  --out DIR  --ne 59 --truth 60 --zmax 16000 --overwrite
Memory: about 31 GB per ensemble held (two at a time), so about 70 GB with the work arrays.
"""
import argparse
import os
import sys
import time
from datetime import datetime, timedelta
from multiprocessing import get_context, shared_memory
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import osse_io as io  # noqa: E402

STORM_THRESH_DBZ = 20.0
TOL = dict(ref=1e-3, qr=1e-7, qs=1e-7, qg=1e-7, T=1e-3, P=1e-1, u=1e-4, v=1e-4, w=1e-4)
DOMAINS = ("global", "reached", "storm", "observed")


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
        """Member files of `exp` valid at t (exp 'prior', 'noda' or an OSSE experiment)."""
        mem = range(1, self.ne + 1)
        if t == self.t0:
            if exp in ("prior", "noda"):
                return [self.gues(m) for m in mem]
            return [self.anal(exp, m) for m in mem]
        return [self.fcst(exp, m, t) for m in mem]

    def truth_file(self, t):
        return self.gues(self.truth) if t == self.t0 else self.fcst("noda", self.truth, t)


# ----------------------------------------------------------------------------------
# Reading an ensemble into shared memory, one process per member file
# ----------------------------------------------------------------------------------
_SHM = {}


def _init(names, shape, box):
    global _SHM
    _SHM = {"box": box, "arr": {}}
    for v, name in names.items():
        shm = shared_memory.SharedMemory(name=name)
        _SHM["arr"][v] = (shm, np.ndarray(shape, np.float32, buffer=shm.buf))


def _read_member(args):
    m, path = args
    x = io.read_wrf(path, _SHM["box"])
    for v, (_, arr) in _SHM["arr"].items():
        arr[m] = x[v]
    return m


class Ensemble:
    """(NE, nk, ny, nx) float32 arrays per variable, in shared memory."""

    def __init__(self, files, box, nproc):
        k0, k1, j0, j1, i0, i1 = box
        shape = (len(files), k1 - k0, j1 - j0, i1 - i0)
        nbytes = int(np.prod(shape)) * 4
        self.shm = {v: shared_memory.SharedMemory(create=True, size=nbytes) for v in io.VARS}
        self.x = {v: np.ndarray(shape, np.float32, buffer=s.buf) for v, s in self.shm.items()}
        missing = [str(f) for f in files if not Path(f).exists()]
        if missing:
            self.close()
            raise FileNotFoundError(f"{len(missing)} member files missing, first {missing[0]}")
        names = {v: s.name for v, s in self.shm.items()}
        with get_context("fork").Pool(nproc, initializer=_init, initargs=(names, shape, box)) as pool:
            for _ in pool.imap_unordered(_read_member, list(enumerate(files))):
                pass

    def close(self):
        self.x = {}
        for s in self.shm.values():
            s.close()
            s.unlink()
        self.shm = {}


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


def sums(st, dom, ref_abserr=None, tol=0.0):
    """Pooled sums of one variable over one domain, total and per level."""
    w = dom.astype(np.float64)
    e = st["err"].astype(np.float64)
    out = dict(n=w.sum(axis=(1, 2)), sse=(w * e * e).sum(axis=(1, 2)), se=(w * e).sum(axis=(1, 2)),
               svar=(w * st["var"]).sum(axis=(1, 2)), scrps=(w * st["crps"]).sum(axis=(1, 2)))
    if ref_abserr is not None:
        s = ref_abserr - np.abs(st["err"])
        out["n_imp"] = ((s > tol) & dom).sum(axis=(1, 2)).astype(np.float64)
        out["n_deg"] = ((s < -tol) & dom).sum(axis=(1, 2)).astype(np.float64)
    return out


# ----------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--osse", default=f"/data/ra000007/{os.environ.get('USER', '')}/osse")
    ap.add_argument("--date", default="20240319")
    ap.add_argument("--exps", default="nt1,nt2,nt4")
    ap.add_argument("--times", default=",".join(f"{18 + m // 60:02d}{m % 60:02d}" for m in range(0, 61, 5)))
    ap.add_argument("--out", default=None)
    ap.add_argument("--ne", type=int, default=59)
    ap.add_argument("--truth", type=int, default=60)
    ap.add_argument("--zmax", type=float, default=16000.0)
    ap.add_argument("--nproc", type=int, default=16)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()

    L = Layout(a.osse, a.date, a.ne, a.truth)
    exps = a.exps.split(",")
    out = Path(a.out) if a.out else L.osse / "osse18" / "verif"
    out.mkdir(parents=True, exist_ok=True)
    times = [datetime.strptime(a.date + hm, "%Y%m%d%H%M") for hm in a.times.split(",")]
    tick = time.time()
    log = lambda msg: print(f"[{time.time() - tick:7.1f} s] {msg}", flush=True)  # noqa: E731

    # ---- verification box and fixed domains -------------------------------------
    mfile = out / "masks.npz"
    if mfile.exists() and not a.overwrite:
        mk = dict(np.load(mfile))
        box = tuple(int(b) for b in mk["box"])
        log(f"masks from {mfile}")
    else:
        if L.t0 not in times:
            sys.exit("masks.npz missing: run 18:00 first")
        info = io.grid_info(L.truth_file(L.t0))
        nlev, nlat, nlon = info["nlev"], info["nlat"], info["nlon"]
        zw = info["zw"][:nlev]
        kmax = max(k for k in range(1, nlev - 2) if zw[k, 2:nlat - 2, 2:nlon - 2].min() <= a.zmax)
        box = (1, kmax + 1, 2, nlat - 2, 2, nlon - 2)
        k0, k1, j0, j1, i0, i1 = box
        region = zw[k0:k1, j0:j1, i0:i1] <= a.zmax
        obs = io.read_radar_obs(L.osse / "osse18" / f"anal_{exps[0]}" / "rad0101.dat")
        kk, jj, ii, _ = io.obs_to_ijk(obs, info)
        inbox = (kk >= k0) & (kk < k1) & (jj >= j0) & (jj < j1) & (ii >= i0) & (ii < i1)
        if not inbox.all():
            log(f"{(~inbox).sum()} observations outside the verification box, not in 'observed'")
        observed = np.zeros(region.shape, bool)
        observed[kk[inbox] - k0, jj[inbox] - j0, ii[inbox] - i0] = True
        zmean = zw[k0:k1, j0:j1, i0:i1].mean(axis=(1, 2))
        mk = dict(box=np.array(box), region=region, observed=observed & region,
                  zmean=zmean, lev5=np.argmin(np.abs(zmean - 5000.0)),
                  lat=info["lat"][j0:j1, i0:i1], lon=info["lon"][j0:j1, i0:i1], n_obs=len(obs))
        log(f"box k {k0}:{k1} j {j0}:{j1} i {i0}:{i1}; {region.sum()} cells, "
            f"{observed.sum()} observed ({len(obs)} observations)")

    for t in times:
        tag = f"{t:%H%M}"
        if (out / f"stats_{tag}.csv").exists() and not a.overwrite:
            log(f"{tag}: exists, skipped")
            continue
        lead = int((t - L.t0).total_seconds() // 60)
        refname = "prior" if t == L.t0 else "noda"
        truth = io.read_wrf(L.truth_file(t), box)
        region = mk["region"]
        colmax = np.where(region, truth["ref"], -np.inf).max(axis=0)
        storm = np.broadcast_to(colmax >= STORM_THRESH_DBZ, region.shape) & region
        log(f"{tag}: truth read, {storm[0].sum()} storm columns")

        ref = Ensemble(L.files(refname, t), box, a.nproc)
        log(f"{tag}: {refname} read")
        if "reached" not in mk:              # 18:00, first experiment against the prior
            first = Ensemble(L.files(exps[0], t), box, a.nproc)
            changed = np.zeros(region.shape, bool)
            for v in io.VARS:
                changed |= (first.x[v] != ref.x[v]).any(axis=0)
            first.close()
            mk["reached"] = changed & region
            j_storm = storm[mk["lev5"]].sum(axis=1)
            mk["jsec"] = int(np.argmax(j_storm))
            np.savez_compressed(mfile, **mk)
            log(f"{tag}: reached {mk['reached'].sum()} cells; masks saved")
        domains = dict(global_=region, reached=mk["reached"], storm=storm, observed=mk["observed"])
        lev5, jsec = int(mk["lev5"]), int(mk["jsec"])

        rows, levels, slices = [], {}, {}
        slices["truth_ref5"] = truth["ref"][lev5]
        slices["truth_colmax"] = np.where(region, truth["ref"], np.nan).max(axis=0)
        slices["truth_sec"] = truth["ref"][:, jsec, :]
        ref_abserr = {}
        for name, ens in [(refname, ref)] + [(e, None) for e in exps]:
            if ens is None:
                try:
                    ens = Ensemble(L.files(name, t), box, a.nproc)
                except FileNotFoundError as err:
                    log(f"{tag}: {name} skipped ({err})")
                    continue
            for v in io.VARS:
                st = field_stats(ens.x[v], truth[v])
                if name == refname:
                    ref_abserr[v] = np.abs(st["err"])
                for dname, dom in domains.items():
                    d = dname.rstrip("_")
                    s = sums(st, dom, None if name == refname else ref_abserr[v], TOL[v])
                    levels[f"{name}|{v}|{d}"] = np.stack([s[k] for k in ("n", "sse", "se", "svar", "scrps")])
                    rows.append(dict(time=tag, lead=lead, exp=name, var=v, domain=d,
                                     **{k: float(np.sum(x)) for k, x in s.items()}))
                if v == "ref":
                    slices[f"{name}_ref5"] = st["mean"][lev5]
                    slices[f"{name}_colmax"] = np.where(region, st["mean"], np.nan).max(axis=0)
                    slices[f"{name}_sec"] = st["mean"][:, jsec, :]
                    slices[f"{name}_spread5"] = np.sqrt(st["var"][lev5])
                    if name != refname:
                        slices[f"{name}_skill5"] = ref_abserr[v][lev5] - np.abs(st["err"][lev5])
                del st
            ens.close()
            log(f"{tag}: {name} done")

        import csv
        with open(out / f"stats_{tag}.csv", "w", newline="") as f:
            cols = ["time", "lead", "exp", "var", "domain", "n", "sse", "se", "svar", "scrps", "n_imp", "n_deg"]
            wr = csv.DictWriter(f, fieldnames=cols)
            wr.writeheader()
            for r in rows:
                wr.writerow({c: r.get(c, "") for c in cols})
        np.savez_compressed(out / f"levels_{tag}.npz", **levels)
        np.savez_compressed(out / f"slices_{tag}.npz", **slices)
        log(f"{tag}: written")


if __name__ == "__main__":
    main()
