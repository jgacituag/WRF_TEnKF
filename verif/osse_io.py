"""Reading WRF output and the synthetic observations for the OSSE verification.

The reflectivity operator is Tong & Xue (2006), method 2 of calc_ref_vr in
src/wrf_common/common_obs_wrf.f90 (the operator the LETKF and synth_obs use), in the
NumPy form of the thesis Chapter 4 notebooks (nbcommon.calc_ref_np). Temperature and
pressure are built as read_var_wrf builds them: P + PB, and (T + 300) * ((P + PB)/P00)
** (Rd/cp) in single precision, with cp = 3.5 Rd.

Python on the Fugaku pre/post nodes: netCDF4 is installed with
    python3.12 -m pip install --no-deps --target /data/ra000007/$USER/pylib netCDF4 cftime certifi
and found through PYLIB (default /data/ra000007/$USER/pylib).
"""
import os
import sys

_PYLIB = os.environ.get("PYLIB", f"/data/ra000007/{os.environ.get('USER', '')}/pylib")
if os.path.isdir(_PYLIB) and _PYLIB not in sys.path:
    sys.path.insert(0, _PYLIB)

import numpy as np               # noqa: E402
from netCDF4 import Dataset      # noqa: E402

PI = 3.141592653589793
RD = 287.0
CP = 3.5 * RD
GG = 9.81
DBZ_MIN = 0.0

# Verified variables. "ref" is reflectivity [dBZ], computed member by member and then
# averaged in dBZ, as in Chapter 4 (mean of H(x), not H of the mean).
VARS = ["ref", "qr", "qs", "qg", "T", "P", "u", "v", "w"]
UNITS = {"ref": "dBZ", "qr": "kg/kg", "qs": "kg/kg", "qg": "kg/kg", "T": "K", "P": "Pa",
         "u": "m/s", "v": "m/s", "w": "m/s"}


def calc_ref(qr, qs, qg, t, p):
    """Equivalent reflectivity Z [mm^6 m^-3], Tong & Xue (2006), method 2 of calc_ref_vr."""
    qr, qs, qg, t, p = (np.asarray(a, np.float64) for a in (qr, qs, qg, t, p))
    nor, nos, nog = 8.0e6, 2.0e6, 4.0e6
    ror, ros, rog, roi = 1000.0, 100.0, 913.0, 917.0
    ki2, kr2 = 0.176, 0.930
    pip = PI ** 1.75
    cf = 1.0e18 * 720 / (pip * nor ** 0.75 * ror ** 1.75)
    cf2 = 1.0e18 * 720 * ki2 * ros ** 0.25 / (pip * kr2 * nos ** 0.75 * roi ** 2)
    cf3 = 1.0e18 * 720 / (pip * nos ** 0.75 * ros ** 1.75)
    cf4 = (1.0e18 * 720 / (pip * nog ** 0.75 * rog ** 1.75)) ** 0.95
    ro = p / (RD * t)
    with np.errstate(invalid="ignore"):
        zr = np.where(qr > 0.0, cf * np.maximum(ro * qr, 0.0) ** 1.75, 0.0)
        zs = np.where(qs > 0.0, np.where(t <= 273.16, cf2, cf3) * np.maximum(ro * qs, 0.0) ** 1.75, 0.0)
        zg = np.where(qg > 0.0, cf4 * np.maximum(ro * qg, 0.0) ** 1.6625, 0.0)
    return zr + zs + zg


def to_dbz(z, dbz_min=DBZ_MIN):
    """dBZ with the LETKF censoring: Z above 10**(dbz_min/10) -> 10 log10 Z, else dbz_min."""
    zmin = 10.0 ** (dbz_min / 10.0)
    out = np.full(np.shape(z), dbz_min, np.float64)
    pos = z > zmin
    out[pos] = 10.0 * np.log10(z[pos])
    return out


def read_wrf(path, box, with_height=False):
    """The verified variables of one wrfout, on mass points, cropped to box.

    box = (k0, k1, j0, j1, i0, i1), Python slices on (bottom_top, south_north, west_east).
    Returns a dict of float32 arrays (k, j, i); with_height adds "z", the height [m] of
    the w level with the same index (PH/g of the cell's lower face), which is the height
    synth_obs gives each observation.
    """
    k0, k1, j0, j1, i0, i1 = box
    s3 = (0, slice(k0, k1), slice(j0, j1), slice(i0, i1))
    with Dataset(path) as nc:
        nc.set_auto_mask(False)
        g = lambda name, sl=s3: np.asarray(nc[name][sl], np.float32)  # noqa: E731
        p = g("P") + g("PB")
        p00 = np.float32(nc["P00"][0]) if "P00" in nc.variables else np.float32(1.0e5)
        if p00 == 0:
            p00 = np.float32(1.0e5)
        t = (g("T") + np.float32(300.0)) * (p / p00) ** np.float32(RD / CP)
        out = {"qr": g("QRAIN"), "qs": g("QSNOW"), "qg": g("QGRAUP"), "T": t, "P": p}
        u = np.asarray(nc["U"][0, k0:k1, j0:j1, i0:i1 + 1], np.float32)
        v = np.asarray(nc["V"][0, k0:k1, j0:j1 + 1, i0:i1], np.float32)
        w = np.asarray(nc["W"][0, k0:k1 + 1, j0:j1, i0:i1], np.float32)
        out["u"] = 0.5 * (u[:, :, :-1] + u[:, :, 1:])
        out["v"] = 0.5 * (v[:, :-1, :] + v[:, 1:, :])
        out["w"] = 0.5 * (w[:-1] + w[1:])
        if with_height:
            out["z"] = ((g("PH") + g("PHB")) / np.float32(GG)).astype(np.float32)
    out["ref"] = to_dbz(calc_ref(out["qr"], out["qs"], out["qg"], out["T"], out["P"])).astype(np.float32)
    return out


def grid_info(path):
    """Grid size, coordinates and P00 of a wrfout."""
    with Dataset(path) as nc:
        nc.set_auto_mask(False)
        info = dict(nlon=nc.dimensions["west_east"].size, nlat=nc.dimensions["south_north"].size,
                    nlev=nc.dimensions["bottom_top"].size,
                    lat=np.asarray(nc["XLAT"][0]), lon=np.asarray(nc["XLONG"][0]),
                    zw=np.asarray((nc["PH"][0] + nc["PHB"][0]) / GG, np.float32))
    return info


def read_radar_obs(path):
    """A radar file written by synth_obs (radar_obs_lonlatz): header lon, lat, z, then
    records (elem, lon, lat, z, value, err, typ). Little- or big-endian, detected from
    the first record marker. Returns a (n, 7) float32 array."""
    raw = open(path, "rb").read()
    endian = "<" if np.frombuffer(raw[:4], "<i4")[0] == 4 else ">"
    i4, f4 = np.dtype(endian + "i4"), np.dtype(endian + "f4")
    pos, recs = 0, []
    nhead = 0
    while pos < len(raw):
        n = int(np.frombuffer(raw[pos:pos + 4], i4)[0])
        body = np.frombuffer(raw[pos + 4:pos + 4 + n], f4)
        pos += 8 + n
        if nhead < 3:
            nhead += 1
            continue
        recs.append(body)
    return np.array(recs, np.float32).reshape(-1, 7)


def obs_to_ijk(obs, info):
    """Grid indices (k, j, i), Python numbering, of observations placed on mass points.

    synth_obs puts each observation at a grid point (i, j) and at z = PH/g of the truth
    at w level k, so the nearest grid point in (lat, lon) and then the level whose
    truth height is nearest give the indices back exactly.
    """
    from scipy.spatial import cKDTree
    lat, lon = info["lat"], info["lon"]
    c = np.cos(np.deg2rad(np.nanmean(lat)))
    tree = cKDTree(np.column_stack([lat.ravel(), lon.ravel() * c]))
    dist, idx = tree.query(np.column_stack([obs[:, 2], obs[:, 1] * c]))
    jj, ii = np.unravel_index(idx, lat.shape)
    kk = np.argmin(np.abs(info["zw"][:, jj, ii] - obs[None, :, 3]), axis=0)
    return kk, jj, ii, dist
