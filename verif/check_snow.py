"""Is the spurious echo of the analyses snow above 0 °C? A diagnostic for 18:00.

For the truth, the prior and each experiment it reads the members one at a time (memory
of one member per reading process) and accumulates, per model level of the verification
box (masks.npz of reduce_osse.py):

  - cells with Z >= THR (35 dBZ by default), split into cells above and below 0 °C,
    averaged over members; the frequency bias per level is that count over the truth's;
  - the mean snow mixing ratio in cells above 0 °C (where the operator treats snow as
    wet snow, with a much larger reflectivity);
  - ensemble means of Z, qs and T along the vertical section of reduce_osse (jsec).

It prints a summary and writes two figures to verif/figures/diag/:
  snow_profiles.png   frequency bias at THR per level, and snow above 0 °C per level;
  snow_section.png    mean qs along the section with the 0 °C line of each ensemble.

Usage (prepost node, after reduce_osse.py --times 1800):
  python3.12 verif/check_snow.py --nproc 8
  options: --osse DIR --verif DIR --exps nt1,nt2,nt4 --ne 59 --truth 60 --thr 35
"""
import argparse
import os
import sys
import time
from multiprocessing import TimeoutError as PoolTimeout, get_context
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import osse_io as io  # noqa: E402
from reduce_osse import Layout, _read_member  # noqa: E402

T0C = 273.15


def accumulate(x, region, thr, jsec, acc):
    """Add one member (dict of (nk, ny, nx) fields) to the accumulators."""
    z, qs, t = x["ref"], x["qs"], x["T"]
    warm = (t > T0C) & region
    ev = (z >= thr) & region
    acc["n_ev"] += ev.sum(axis=(1, 2))
    acc["n_ev_warm"] += (ev & warm).sum(axis=(1, 2))
    acc["qs_warm"] += np.where(warm, qs, 0).sum(axis=(1, 2), dtype=np.float64)
    acc["n_warm"] += warm.sum(axis=(1, 2))
    for v in ("ref", "qs", "T"):
        acc[f"sec_{v}"] += x[v][:, jsec, :]
    acc["n"] += 1


def new_acc(nk, nx):
    a = {k: np.zeros(nk) for k in ("n_ev", "n_ev_warm", "qs_warm", "n_warm")}
    a.update({f"sec_{v}": np.zeros((nk, nx)) for v in ("ref", "qs", "T")})
    a["n"] = 0
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--osse", default=f"/data/ra000007/{os.environ.get('USER', '')}/osse")
    ap.add_argument("--verif", default=None, help="folder of masks.npz [<osse>/osse18/verif]")
    ap.add_argument("--date", default="20240319")
    ap.add_argument("--exps", default="nt1,nt2,nt4")
    ap.add_argument("--ne", type=int, default=59)
    ap.add_argument("--truth", type=int, default=60)
    ap.add_argument("--thr", type=float, default=35.0)
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--timeout", type=float, default=900.0)
    ap.add_argument("--figdir", default=str(Path(__file__).resolve().parent / "figures" / "diag"))
    a = ap.parse_args()

    L = Layout(a.osse, a.date, a.ne, a.truth)
    verif = Path(a.verif) if a.verif else L.osse / "osse18" / "verif"
    mk = dict(np.load(verif / "masks.npz"))
    box = tuple(int(b) for b in mk["box"])
    region, jsec, zkm = mk["region"], int(mk["jsec"]), mk["zmean"] / 1e3
    nk, nx = region.shape[0], region.shape[2]
    tick = time.time()
    log = lambda msg: print(f"[{time.time() - tick:7.1f} s] {msg}", flush=True)  # noqa: E731
    nproc = min(a.nproc, len(os.sched_getaffinity(0)))
    pool = get_context("spawn").Pool(nproc)

    res = {}
    try:
        tr = new_acc(nk, nx)
        accumulate(io.read_wrf(L.truth_file(L.t0), box), region, a.thr, jsec, tr)
        res["truth"] = tr
        log("truth read")
        for name in ["prior"] + a.exps.split(","):
            files = L.files(name, L.t0)
            if not all(Path(f).exists() for f in files):
                log(f"{name}: member files missing, skipped")
                continue
            acc = new_acc(nk, nx)
            it = pool.imap_unordered(_read_member, [(m, f, box) for m, f in enumerate(files)])
            for _ in files:
                try:
                    _, x = it.next(timeout=a.timeout)
                except PoolTimeout as err:
                    raise RuntimeError(f"{name}: no member in {a.timeout} s") from err
                accumulate(x, region, a.thr, jsec, acc)
                del x
            res[name] = acc
            log(f"{name}: {acc['n']} members")
    finally:
        pool.terminate()

    # ---- summary -------------------------------------------------------------------
    t = res["truth"]
    ncell = region.sum(axis=(1, 2))
    warm_frac = t["n_warm"] / np.maximum(ncell, 1)
    k0c = int(np.argmin(np.abs(warm_frac - 0.5)))
    print(f"\n0 °C level, domain average (it varies with latitude; the counts use each cell's own T): about {zkm[k0c]:.2f} km (level {k0c})")
    print(f"Cells with Z >= {a.thr:g} dBZ, total over levels (ensemble mean count per member):")
    print(f"{'':8s}{'total':>12s}{'above 0 °C':>12s}{'fraction':>10s}{'bias':>8s}"
          f"{'snow above 0 °C [g/kg]':>24s}")
    for name, r in res.items():
        n = r["n_ev"].sum() / r["n"]
        nw = r["n_ev_warm"].sum() / r["n"]
        qsw = 1e3 * r["qs_warm"].sum() / max(r["n_warm"].sum(), 1)
        print(f"{name:8s}{n:12.0f}{nw:12.0f}{nw / max(n, 1):10.2f}"
              f"{n / max(t['n_ev'].sum(), 1):8.2f}{qsw:24.4f}")

    # ---- figures -------------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import vcommon as vc
    figdir = Path(a.figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    col = lambda n: "dimgray" if n == "truth" else vc.exp_color(n)  # noqa: E731
    lab = lambda n: "truth" if n == "truth" else vc.exp_label(n, "en")  # noqa: E731

    fig, axs = plt.subplots(1, 2, figsize=(vc.PAGE_W * 0.8, 0.45 * vc.PAGE_W), sharey=True)
    for name, r in res.items():
        with np.errstate(invalid="ignore", divide="ignore"):
            bias = (r["n_ev"] / r["n"]) / t["n_ev"]
            qsw = 1e3 * r["qs_warm"] / r["n_warm"]
        if name != "truth":
            axs[0].plot(bias, zkm, color=col(name), label=lab(name))
        axs[1].plot(qsw, zkm, color=col(name), label=lab(name),
                    ls="--" if name == "truth" else "-")
    for ax in axs:
        ax.axhline(zkm[k0c], color="k", lw=0.6, ls=":")
    axs[0].axvline(1, color="k", lw=0.6)
    axs[0].set_xscale("log")
    axs[0].set_xlabel(f"frequency bias, Z ≥ {a.thr:g} dBZ")
    axs[1].set_xlabel("mean qs above 0 °C [g/kg]")
    axs[0].set_ylabel("height [km]")
    axs[0].set_ylim(0, 8)
    axs[1].legend(loc="best")
    vc.panel_label(axs[0], "a"); vc.panel_label(axs[1], "b")
    fig.tight_layout()
    fig.savefig(figdir / "snow_profiles.png", dpi=200)
    plt.close(fig)

    xlon = mk["lon"][jsec, :]
    names = list(res)
    fig, axs = plt.subplots(len(names), 1, figsize=(vc.PAGE_W * 0.7, 0.2 * vc.PAGE_W * len(names)),
                            sharex=True)
    for ax, name in zip(axs, names):
        r = res[name]
        qs = 1e3 * r["sec_qs"] / r["n"]
        pm = ax.pcolormesh(xlon, zkm, qs, cmap="Blues", vmin=0, vmax=2, shading="auto", rasterized=True)
        ax.contour(xlon, zkm, r["sec_T"] / r["n"], levels=[T0C], colors="r", linewidths=0.8)
        ax.contour(xlon, zkm, r["sec_ref"] / r["n"], levels=[a.thr], colors="k", linewidths=0.5)
        ax.set_ylim(0, 8)
        ax.set_ylabel("height [km]")
        ax.text(0.01, 0.95, lab(name), transform=ax.transAxes, va="top", fontsize=8,
                bbox=dict(fc="white", ec="none", alpha=0.7))
    axs[-1].set_xlabel("longitude")
    fig.colorbar(pm, ax=axs, shrink=0.6, label="mean qs [g/kg]")
    fig.suptitle(f"latitude {np.nanmean(mk['lat'][jsec]):.2f}°; red 0 °C, black {a.thr:g} dBZ "
                 f"(ensemble means)", y=1.0, fontsize=8)
    fig.savefig(figdir / "snow_section.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    log(f"figures in {figdir}")


if __name__ == "__main__":
    main()
