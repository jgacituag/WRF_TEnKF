"""Shared helpers for the OSSE verification notebooks (V1, ...).

Style, palettes, vocabulary, number formatting and bilingual saving follow the
Chapter 4 notebooks (WRF_Single_Cycle_Assimilation/Notebooks/nbcommon.py): the same
rcParams, dBZ scale, diverging skill map, Nt shades and save(draw, stem) convention,
copied here so that this repository does not depend on that one.

Metric definitions, as in Chapter 4 (da/metrics.py):
  RMSE   sqrt(sum err^2 / n), err = ensemble mean - truth (reflectivity: mean of H(x_i))
  spread sqrt((Ne+1)/Ne * sum var / n), var with ddof = 1
  CRPS   sum CRPS / n
  SS     1 - M_exp / M_ref, with M pooled over the domain before the ratio
  F      fraction of the cells that change beyond TOL whose |error| decreases
Spatial scores on column-maximum reflectivity (spatial_scores), event Z >= threshold:
  FSS    fractions skill score (Roberts & Lean 2008), 1 - sum (P - O)^2 / (sum P^2 + sum O^2),
         with O the fraction of truth columns with the event in a square of side s and P
         the same fraction averaged over the members (the ensemble fraction)
  AUC    area under the ROC curve of the neighborhood ensemble probability (the fraction
         of members with the event within NEP_R grid points), against the truth event at
         the column, as Snook et al. (2015) and Supinie et al. (2016)
  bias   frequency bias, member event area over truth event area, averaged over members
"""
import os
import re as _re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap

HERE = Path(__file__).resolve().parent
VERIF = Path(os.environ.get("VERIF_DIR", f"/data/ra000007/{os.environ.get('USER', '')}/osse/osse18/verif"))
FIG_DIR = Path(os.environ.get("VERIF_FIG_DIR", HERE / "figures"))
NE = 59

# ---------------------------------------------------------------------------------
# Style (nbcommon.set_style)
# ---------------------------------------------------------------------------------
MM = 1 / 25.4
PAGE_W = 174 * MM


def set_style():
    plt.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 600, "savefig.bbox": "tight",
        "font.size": 9, "font.family": "serif", "mathtext.fontset": "stix",
        "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "axes.linewidth": 0.6,
        "lines.linewidth": 1.4, "axes.grid": True, "grid.linewidth": 0.4, "grid.alpha": 0.30,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "figure.facecolor": "white", "axes.facecolor": "white",
    })


set_style()


def panel_label(ax, letter, loc="upper left", pad=0.02, fontsize=9):
    x, ha = (pad, "left") if "left" in loc else (1 - pad, "right")
    y, va = (1 - pad, "top") if "upper" in loc else (pad, "bottom")
    ax.text(x, y, f"({letter})", transform=ax.transAxes, ha=ha, va=va, fontweight="bold",
            fontsize=fontsize, bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))


DBZ_LEVELS = np.arange(0, 66, 5)
DBZ_CMAP = "Spectral_r"
DBZ_NORM = BoundaryNorm(DBZ_LEVELS, plt.get_cmap(DBZ_CMAP).N)
CMAP_DIV = "RdBu_r"           # red = positive = better, as CMAP_SKILL in nbcommon
_NT_SHADES = ["#9ec5f4", "#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"]


def exp_color(exp):
    if exp in ("prior", "noda"):
        return "black"
    return _NT_SHADES[min(int(exp[2:]), 5) - 1]


def exp_label(exp, lang):
    if exp == "prior":
        return "prior"
    if exp == "noda":
        return {"es": "sin asimilación", "en": "no DA"}[lang]
    return f"T{int(exp[2:])}"


def nt_of(exp):
    return int(exp[2:])


# ---------------------------------------------------------------------------------
# Vocabulary (nbcommon.VOCAB, plus the forecast terms)
# ---------------------------------------------------------------------------------
LANGS = ("es", "en")
SHOW = "es"
FIG_DPI = 300
VOCAB = {
    "es": {
        "ss_rmse": r"$\mathrm{SS}_{\mathrm{RMSE}}$", "ss_crps": r"$\mathrm{SS}_{\mathrm{CRPS}}$",
        "frac_imp": "fracción que mejora", "nt": r"pasos de temperado $N_t$",
        "truth": "verdad", "ens_mean": "media del ensamble", "ens_spread": "dispersión",
        "dom_global": "dominio", "dom_reached": "celdas alcanzadas", "dom_storm": "tormenta",
        "dom_observed": "celdas observadas",
        "dom_echo": "celdas con eco", "fss": "FSS", "auc": "AUC (ROC)",
        "freq_bias": "sesgo de frecuencia", "scale_km": "escala [km]", "thr": "umbral", "lead": "plazo [min]", "height": "altura [km]",
        "rmse": "RMSE", "spread_ratio": "dispersión / RMSE", "ref": "reflectividad [dBZ]",
        "skill": r"$|e_{\mathrm{ref}}| - |e|$ [dBZ]", "x_km": "x [km]",
        "lon": "longitud", "lat": "latitud",
    },
    "en": {
        "ss_rmse": r"$\mathrm{SS}_{\mathrm{RMSE}}$", "ss_crps": r"$\mathrm{SS}_{\mathrm{CRPS}}$",
        "frac_imp": "fraction improving", "nt": r"tempering steps $N_t$",
        "truth": "truth", "ens_mean": "ensemble mean", "ens_spread": "spread",
        "dom_global": "domain", "dom_reached": "reached cells", "dom_storm": "storm",
        "dom_observed": "observed cells",
        "dom_echo": "echo cells", "fss": "FSS", "auc": "ROC AUC",
        "freq_bias": "frequency bias", "scale_km": "scale [km]", "thr": "threshold", "lead": "lead time [min]", "height": "height [km]",
        "rmse": "RMSE", "spread_ratio": "spread / RMSE", "ref": "reflectivity [dBZ]",
        "skill": r"$|e_{\mathrm{ref}}| - |e|$ [dBZ]", "x_km": "x [km]",
        "lon": "longitude", "lat": "latitude",
    },
}
assert set(VOCAB["es"]) == set(VOCAB["en"])
DOMAINS = ("global", "echo", "observed")

_DECIMAL_POINT = _re.compile(r"(?<=\d)\.(?=\d)")
_WEST = _re.compile("(°|\\^\\{\\\\circ\\}\\$?)W\\b")


def _check_lang(lang):
    if lang not in LANGS:
        raise ValueError(f"language {lang!r}; expected one of {LANGS}")


def _to_spanish(s):
    return _WEST.sub(r"\1O", _DECIMAL_POINT.sub(",", s))


def fmt(x, lang, spec="g"):
    """A number as text in language `lang` (decimal comma in Spanish)."""
    _check_lang(lang)
    s = format(x, spec)
    if plt.rcParams["axes.unicode_minus"]:
        s = s.replace("-", "−")
    if lang == "es":
        s = s.replace(",", "\x00").replace(".", ",").replace("\x00", " ")
    return s


class _SpanishFormatter(mticker.Formatter):
    def __init__(self, base):
        self.base = base

    def set_axis(self, axis):
        super().set_axis(axis)
        self.base.set_axis(axis)

    def set_locs(self, locs):
        self.base.set_locs(locs)

    def __call__(self, x, pos=None):
        return _to_spanish(self.base(x, pos))

    def format_ticks(self, values):
        return [_to_spanish(s) for s in self.base.format_ticks(values)]

    def get_offset(self):
        return _to_spanish(self.base.get_offset())


def localize(fig, lang):
    _check_lang(lang)

    def swap(f):
        wrapped = isinstance(f, _SpanishFormatter)
        if lang == "es" and not wrapped:
            return _SpanishFormatter(f)
        if lang == "en" and wrapped:
            return f.base
        return f

    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            axis.set_major_formatter(swap(axis.get_major_formatter()))
            axis.set_minor_formatter(swap(axis.get_minor_formatter()))
    return fig


def save(fig, stem, show=None, dpi=FIG_DPI, langs=None):
    """Draw and save a figure in both languages: fig is a function lang -> Figure.
    Writes FIG_DIR/<lang>/<stem>.png and leaves the `show` language open (nbcommon.save)."""
    if not callable(fig) or isinstance(fig, plt.Figure):
        raise TypeError("save() expects a function lang -> Figure")
    show = SHOW if show is None else show
    langs = LANGS if langs is None else tuple(langs)
    paths = {}
    for lang in langs:
        f = fig(lang)
        localize(f, lang)
        p = FIG_DIR / lang / f"{stem}.png"
        p.parent.mkdir(parents=True, exist_ok=True)
        f.savefig(p, dpi=dpi, bbox_inches="tight")
        paths[lang] = p
        if lang != show:
            plt.close(f)
    print("  ".join(str(p) for p in paths.values()))
    return paths


# ---------------------------------------------------------------------------------
# Loading the reduced files (verif/reduce_osse.py)
# ---------------------------------------------------------------------------------
def load_stats(verif=VERIF):
    """All stats_<HHMM>.csv as one DataFrame, with the derived metrics."""
    files = sorted(Path(verif).glob("stats_*.csv"))
    if not files:
        raise FileNotFoundError(f"no stats_*.csv in {verif}: run verif/reduce_osse.py")
    df = pd.concat([pd.read_csv(f, dtype={"time": str}) for f in files], ignore_index=True)
    df["rmse"] = np.sqrt(df.sse / df.n)
    df["bias"] = df.se / df.n
    df["spread"] = np.sqrt((NE + 1) / NE * df.svar / df.n)
    df["crps"] = df.scrps / df.n
    df["frac_imp"] = df.n_imp / (df.n_imp + df.n_deg)
    return df


def skill(df):
    """Skill scores of every experiment against the reference of its time (prior at
    18:00, the forecast without DA after), pooled before the ratio."""
    ref = df[df.exp.isin(["prior", "noda"])].set_index(["time", "var", "domain"])
    exp = df[~df.exp.isin(["prior", "noda"])].copy()
    key = list(zip(exp.time, exp["var"], exp.domain))
    with np.errstate(invalid="ignore", divide="ignore"):     # variables absent from the truth
        exp["ss_rmse"] = 1 - exp.rmse.values / ref.loc[key, "rmse"].values
        exp["ss_crps"] = 1 - exp.crps.values / ref.loc[key, "crps"].values
    exp["rmse_ref"] = ref.loc[key, "rmse"].values
    exp["crps_ref"] = ref.loc[key, "crps"].values
    return exp


def load_npz(name, verif=VERIF):
    with np.load(Path(verif) / name) as z:
        return {k: z[k] for k in z.files}


def experiments(df):
    return sorted({e for e in df.exp.unique() if e not in ("prior", "noda")}, key=nt_of)


# ---------------------------------------------------------------------------------
# Spatial scores of column-maximum reflectivity (colmax_<HHMM>.npz)
# ---------------------------------------------------------------------------------
THRESHOLDS = (20.0, 35.0)       # dBZ
SCALES = (1, 3, 7, 15)          # FSS square side in grid points (4, 12, 28, 60 km at 4 km)
NEP_R = 1                       # neighborhood radius of the ensemble probability [grid points]
DX_KM = 4.0


def _disk(r):
    y, x = np.mgrid[-r:r + 1, -r:r + 1]
    return x * x + y * y <= r * r


def fss(p, o, s):
    """FSS of forecast and observed fields of 0/1 or fractions, square side s."""
    from scipy.ndimage import uniform_filter
    pf = uniform_filter(p.astype(np.float64), s, mode="constant")
    of = uniform_filter(o.astype(np.float64), s, mode="constant")
    den = (pf ** 2).sum() + (of ** 2).sum()
    return np.nan if den == 0 else 1 - ((pf - of) ** 2).sum() / den


def roc_auc(prob, event):
    """Area under the ROC curve; prob takes the values k/NE."""
    prob, event = prob.ravel(), event.ravel()
    ne_ev, ne_no = event.sum(), (~event).sum()
    if ne_ev == 0 or ne_no == 0:
        return np.nan
    levels = np.unique(prob)[::-1]
    pod = [0.0] + [(event & (prob >= q)).sum() / ne_ev for q in levels]
    pofd = [0.0] + [(~event & (prob >= q)).sum() / ne_no for q in levels]
    pod.append(1.0); pofd.append(1.0)
    return float(np.trapezoid(pod, pofd) if hasattr(np, "trapezoid") else np.trapz(pod, pofd))


def spatial_scores(verif=VERIF, thresholds=THRESHOLDS, scales=SCALES, nep_r=NEP_R):
    """FSS, ROC AUC and frequency bias of every ensemble in every colmax_<HHMM>.npz, as a
    long DataFrame (time, lead, exp, thr, metric, scale_km, value)."""
    from scipy.ndimage import maximum_filter
    rows = []
    files = sorted(Path(verif).glob("colmax_*.npz"))
    if not files:
        raise FileNotFoundError(f"no colmax_*.npz in {verif}: run verif/reduce_osse.py")
    for f in files:
        tag = f.stem.split("_")[1]
        lead = (int(tag[:2]) - 18) * 60 + int(tag[2:])
        with np.load(f) as z:
            truth = z["truth"].astype(np.float32)
            for exp in [k for k in z.files if k != "truth"]:
                ens = z[exp].astype(np.float32)
                for thr in thresholds:
                    obs = truth >= thr
                    ev = ens >= thr
                    base = dict(time=tag, lead=lead, exp=exp, thr=thr)
                    p = ev.mean(axis=0)
                    for s in scales:
                        rows.append(dict(base, metric="fss", scale_km=s * DX_KM, value=fss(p, obs, s)))
                    nb = np.stack([maximum_filter(m, footprint=_disk(nep_r), mode="constant")
                                   for m in ev]).mean(axis=0)
                    rows.append(dict(base, metric="auc", scale_km=nep_r * DX_KM, value=roc_auc(nb, obs)))
                    nobs = obs.sum()
                    rows.append(dict(base, metric="bias", scale_km=0.0,
                                     value=ev.sum(axis=(1, 2)).mean() / nobs if nobs else np.nan))
    return pd.DataFrame(rows)
