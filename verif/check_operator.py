"""Check the NumPy reflectivity operator against the Fortran one.

synth_obs writes, beside the observations, the truth's own reflectivity at every
observation point (truth_rad0101.dat), computed with calc_ref_vr in Fortran. This
script recomputes it with verif/osse_io.py at the same grid points and reports the
differences in dBZ.

Usage: python3.12 verif/check_operator.py TRUTH_WRFOUT TRUTH_RAD_FILE
  e.g. .../osse/da4km/HIST/GUES/20240319180000/gues00060 .../osse/osse18/anal_nt1/truth_rad0101.dat
"""
import sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import osse_io as io  # noqa: E402

truth_file, rad_file = sys.argv[1], sys.argv[2]
info = io.grid_info(truth_file)
obs = io.read_radar_obs(rad_file)
k, j, i, dist = io.obs_to_ijk(obs, info)
box = (0, info["nlev"], 0, info["nlat"], 0, info["nlon"])
x = io.read_wrf(truth_file, box)
mine = x["ref"][k, j, i].astype(np.float64)
fortran = io.to_dbz(obs[:, 4].astype(np.float64))
d = mine - fortran
echo = fortran > io.DBZ_MIN
print(f"{len(obs)} observations; grid match distance max {dist.max():.2e} (deg)")
print(f"truth echo at {echo.sum()} of them")
print(f"|numpy - fortran| [dBZ]: max {np.abs(d).max():.3e}, mean {np.abs(d).mean():.3e}, "
      f"99.9 % {np.quantile(np.abs(d), 0.999):.3e}")
print(f"echo/no-echo disagreements: {np.sum((mine > io.DBZ_MIN) != echo)}")
ok = np.abs(d).max() < 0.01
print("PASSED" if ok else "CHECK: differences above 0.01 dBZ")
sys.exit(0 if ok else 1)
