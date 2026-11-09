"""Write small synthetic wrfout files for local tests of the LETKF tools.

Not physics: just enough variables, attributes and structure for
set_common_wrf / read_grd, on a Mercator grid, with a rain blob whose
position and intensity change with the member. Needs numpy and netCDF4
(local machine only; nothing here runs on the cluster).

Usage: python3 fake_wrfout.py OUTDIR NMEM
Writes OUTDIR/mem00001 ... OUTDIR/mem<NMEM>.
"""
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

NX, NY, NZ, NSOIL = 40, 30, 20, 4          # mass points
DX = 2000.0
LAT1, LON1, TRUELAT, STDLON = -36.0, -60.0, -35.0, -60.0
RE = 6371.3e3                              # same radius as common.f90


def mercator_latlon(i, j):
    """WRF Mercator (map_utils set_merc / ij_to_latlon), 1-based i, j."""
    dlon = DX / (RE * np.cos(np.deg2rad(TRUELAT)))
    rsw = np.log(np.tan(0.5 * np.deg2rad(LAT1 + 90.0))) / dlon
    lat = 2.0 * np.rad2deg(np.arctan(np.exp(dlon * (rsw + j - 1.0)))) - 90.0
    lon = LON1 + np.rad2deg((i - 1.0) * dlon)
    return lat, lon


def write_member(path, imem, rng):
    nc = Dataset(path, "w", format="NETCDF4")
    for name, size in (("Time", None), ("DateStrLen", 19), ("west_east", NX),
                       ("west_east_stag", NX + 1), ("south_north", NY),
                       ("south_north_stag", NY + 1), ("bottom_top", NZ),
                       ("bottom_top_stag", NZ + 1), ("soil_layers_stag", NSOIL)):
        nc.createDimension(name, size)
    nc.setncatts({"DX": DX, "DY": DX, "MAP_PROJ": np.int32(3), "TRUELAT1": TRUELAT,
                  "TRUELAT2": TRUELAT, "STAND_LON": STDLON})

    def var(name, dims, data):
        v = nc.createVariable(name, "f4", dims)
        v[:] = data
        return v

    times = nc.createVariable("Times", "S1", ("Time", "DateStrLen"))
    times[0, :] = np.array(list("2024-03-19_18:00:00"), dtype="S1")

    m3 = ("Time", "bottom_top", "south_north", "west_east")
    jj, ii = np.meshgrid(np.arange(1, NY + 1), np.arange(1, NX + 1), indexing="ij")
    lat, lon = mercator_latlon(ii.astype(float), jj.astype(float))
    ju, iu = np.meshgrid(np.arange(1, NY + 1.0), np.arange(0.5, NX + 1.0), indexing="ij")
    latu, lonu = mercator_latlon(iu, ju)
    jv, iv = np.meshgrid(np.arange(0.5, NY + 1.0), np.arange(1, NX + 1.0), indexing="ij")
    latv, lonv = mercator_latlon(iv, jv)

    zw = 250.0 * np.arange(NZ + 1)                       # w levels [m]
    zm = 0.5 * (zw[1:] + zw[:-1])
    pb = 1.0e5 * np.exp(-zm / 8000.0)

    # Rain blob: centre and amplitude change with the member
    ci = NX / 2 + rng.normal(0, 3)
    cj = NY / 2 + rng.normal(0, 3)
    amp = 2.0e-3 * rng.uniform(0.2, 1.0) if imem % 4 else 0.0   # some members dry
    r2 = ((ii - ci) ** 2 + (jj - cj) ** 2) / 4.0 ** 2
    blob2d = amp * np.exp(-r2)
    vert = np.exp(-((zm - 2000.0) / 1200.0) ** 2)
    qr = blob2d[None, :, :] * vert[:, None, None]

    var("U", ("Time", "bottom_top", "south_north", "west_east_stag"), np.zeros((1, NZ, NY, NX + 1)))
    var("V", ("Time", "bottom_top", "south_north_stag", "west_east"), np.zeros((1, NZ, NY + 1, NX)))
    var("W", ("Time", "bottom_top_stag", "south_north", "west_east"), np.zeros((1, NZ + 1, NY, NX)))
    # Potential temperature perturbation: a smooth blob of random amplitude,
    # so that temperature observations see ensemble spread.
    tamp = rng.normal(0.0, 1.0)
    rt2 = ((ii - NX / 3) ** 2 + (jj - NY / 3) ** 2) / 8.0 ** 2
    tvert = np.exp(-((zm - 1500.0) / 1500.0) ** 2)
    var("T", m3, (tamp * np.exp(-rt2))[None, None, :, :] * tvert[None, :, None, None])
    var("P", m3, np.zeros((1, NZ, NY, NX)))
    var("PB", m3, np.broadcast_to(pb[None, :, None, None], (1, NZ, NY, NX)))
    w3 = ("Time", "bottom_top_stag", "south_north", "west_east")
    var("PH", w3, np.zeros((1, NZ + 1, NY, NX)))
    var("PHB", w3, np.broadcast_to(9.81 * zw[None, :, None, None], (1, NZ + 1, NY, NX)))
    var("QVAPOR", m3, np.full((1, NZ, NY, NX), 0.01))
    for name in ("QCLOUD", "QICE", "QSNOW", "QGRAUP"):
        var(name, m3, np.zeros((1, NZ, NY, NX)))
    var("QRAIN", m3, qr[None])
    m2 = ("Time", "south_north", "west_east")
    for name, val in (("PSFC", 1.0e5), ("T2", 300.0), ("Q2", 0.01), ("U10", 0.0), ("V10", 0.0),
                      ("LANDMASK", 1.0), ("TSK", 300.0), ("SNOW", 0.0), ("SNOWH", 0.0),
                      ("SEAICE", 0.0), ("CANWAT", 0.0), ("HGT", 0.0), ("COSALPHA", 1.0),
                      ("SINALPHA", 0.0), ("MU", 0.0), ("MUB", 9.0e4)):
        var(name, m2, np.full((1, NY, NX), val))
    var("XLAT", m2, lat[None]); var("XLONG", m2, lon[None])
    var("XLAT_U", ("Time", "south_north", "west_east_stag"), latu[None])
    var("XLONG_U", ("Time", "south_north", "west_east_stag"), lonu[None])
    var("XLAT_V", ("Time", "south_north_stag", "west_east"), latv[None])
    var("XLONG_V", ("Time", "south_north_stag", "west_east"), lonv[None])
    s3 = ("Time", "soil_layers_stag", "south_north", "west_east")
    var("SMOIS", s3, np.full((1, NSOIL, NY, NX), 0.3))
    var("TSLB", s3, np.full((1, NSOIL, NY, NX), 290.0))
    var("P_TOP", ("Time",), [1.0e5 * np.exp(-zw[-1] / 8000.0)])   # pressure of the top w level
    var("P00", ("Time",), [1.0e5])
    var("DNW", ("Time", "bottom_top"), np.full((1, NZ), -1.0 / NZ))
    var("ZNU", ("Time", "bottom_top"), np.linspace(1, 0, NZ)[None])
    nc.close()


def main():
    out, nmem = Path(sys.argv[1]), int(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(12345)
    for m in range(1, nmem + 1):
        write_member(out / f"mem{m:05d}", m, rng)


if __name__ == "__main__":
    main()
