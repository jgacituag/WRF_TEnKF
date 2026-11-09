# WRF_TEnKF

Tempered LETKF (TEnKF, LETKF-T) for WRF, built on the WRF-LETKF code of
Takemasa Miyoshi and Juan Ruiz.

## Origin and authorship

This repository is a derivative of the WRF-LETKF system distributed in
[`wrfutil_hydra`](https://github.com/jgacituag/wrfutil_hydra), commit
`586cc96`. Nearly all of the code here was written by others:

- **Takemasa Miyoshi** wrote the LETKF code (2009 onward), based on
  Ott et al. (2004) and Hunt et al. (2007), including the MPI framework
  and direct NetCDF I/O.
- **Juan Ruiz** adapted it to WRF and extended it: radar forward
  operators, parameter estimation and the namelist interface
  (`common_namelist.f90`). He maintains `wrfutil_hydra`.
- Changes marked `!GYL` in the source (RTPP/RTPS relaxation through the
  LETKF weights, minimum inflation) are by Guo-Yuan Lien.
- `run_LETKF.sh` comes from the Servicio Meteorológico Nacional
  (Maximiliano A. Sacco, Yanina García Skabar, María Eugenia Dillon,
  Cynthia Mariana Matsudo).
- `SFMT.f90` is Atsushi Ito's Fortran port of SFMT by Mutsuo Saito and
  Makoto Matsumoto. Other routines name their authors in their headers.

Jorge Gacitúa adds likelihood tempering on top of this code. Every change
with respect to the original is visible with

```
git diff reference
```

where the tag `reference` holds the files exactly as they were in
`wrfutil_hydra@586cc96` (only the directory layout and the build script
are new).

The code is distributed under the GNU General Public License v3, the
license of the original (see `LICENSE`).

## Method

The analysis is split into `ntemp` steps with weights α_i that sum to one.
Step i assimilates the same observations with error variance R/α_i, and
recomputes H(x) on the ensemble updated by the previous step. With
`ntemp = 1` the code must reproduce the original LETKF bit for bit; the
regression test in `tests/` checks this and must pass before any
experiment.

Weights follow the exponential scheme w_i ∝ exp(−α_s (Nt+1)/i),
normalized to sum to one. α_s = 0 gives equal weights; larger α_s moves
the weight to the last steps.

## Namelist

```
&tempering
  ntemp   = 1        ! 1 = original LETKF
  alpha_s = 1.0d0    ! delay parameter of the exponential weights
  alpha_list = ...   ! optional explicit weights (override alpha_s)
/
```

A missing `&tempering` block leaves `ntemp = 1`.

## Layout

```
src/common/      generic LETKF modules (Miyoshi)
src/wrf_common/  WRF interface, MPI I/O, observation operators, namelist
src/letkf/       main program, analysis (letkf_tools), observations (letkf_obs)
namelists/       letkf.namelist template used by bin/run_LETKF.sh
bin/             run_LETKF.sh for wrfutil, with the tempering placeholders
tests/           ntemp = 1 regression test against the original executable
make.sh          build script; writes dist/letkf_<COMPILER>.tar
```

## Tests

- `tests/test_weights.sh [COMPILER]` runs serially, without WRF data. It
  checks the weights against thesis Table 2.1 and the reading of the
  `&tempering` block (absent, present, misspelled, weights not summing
  to one).
- `tests/run_regression.sh CASE_DIR REF_EXE NEW_EXE WORK_DIR` runs on
  Fugaku, inside one job. It runs the original executable twice and this
  one once with `ntemp = 1`, on fresh copies of the same inputs, and
  compares every file `letkf.exe` writes with `cmp`. Exit code 0 means
  bit-for-bit identical, 1 means different, 2 means the original code is
  not reproducible in that layout. Run it after every change to `src/`
  and before every experiment.

- `tests/make_case.sh EXPDIR DATE NBV CASE_DIR` builds a test case from a
  saved wrfutil experiment, with links only. wrfutil keeps no prior
  members, so the case uses the saved analyses as "prior": meaningless
  physically, identical for every executable.
- `tests/jobs/regression_fugaku.sh` is the PJM job for Fugaku. It
  compares the `reference` build against `main` (the source change) and,
  once, the `letkf.exe` of the experiments against the `reference` build
  (that our build is equivalent). Instructions are in its header.

Both executables of a comparison must come from the same `make.sh` and
compiler flags; otherwise a difference may come from the build and not
from the source.

## OSSE tools

- `radar_obs_lonlatz` (`&radar_da`, false by default): radar records hold
  longitude, latitude and height above sea level [m] instead of azimuth,
  elevation and range, and `aer2llz` is skipped. Reflectivity only.
- `synth_obs.exe` (`src/osse/synth_obs.f90`) writes synthetic reflectivity
  on model grid points from a truth member, in that format. It uses the
  analysis operator (`Trans_XtoY`), keeps the points where at least
  `min_active` ensemble members have echo, and adds N(0, `obs_err`²) noise
  with a fixed seed. Settings in `synth_obs.namelist`, documented in the
  source header.
- `namelists/letkf_osse_cap4.namelist`: LETKF settings of the thesis
  Chapter 4 multi-obs experiments (2 km localization, no inflation, no
  RTPP, no reflectivity gross-error check, network of points with at
  least one member with echo).
- `osse/make_fcst_exp.sh` builds a wrfutil DAFCST experiment that
  forecasts from a given set of analyses, reusing the scripts, namelists,
  WRF executables and met_em of an existing experiment (copied or linked,
  never modified). `--dry-run` checks every input first.
- `osse/make_da_exp.sh` builds a wrfutil DA-cycle experiment that repeats
  an existing one over a new period (same scripts, namelists, executables,
  GEFS boundaries and radar data; only the dates, `BASEDIR` and
  `EXPTYPE=DACYCLE` change). `--dry-run` checks the executables, the GEFS
  files of every step and the radar files; `--smoke` makes a 4-member,
  8-node version for a check in the interactive group.
- `tests/local/test_osse_tools.sh` runs the chain on small synthetic
  wrfout files on a workstation (needs numpy and netCDF4; nothing of it
  runs on the cluster).

## Use with wrfutil

On Fugaku, keep the repository, builds, test cases and results in the
group data area (`/data/ra000007/$USER`, on `/vol0004`), not in the home
area (20 GiB per user): jobs must be submitted from the data area and read
their files from it. Directories that will receive files larger than 1 GB
from compute nodes (wrfout) need striping, `lfs setstripe -c 8 <dir>`.

`wrfutil_hydra` stays unchanged. In the experiment configuration, point
`LETKFDIR` to a directory holding `dist/letkf_<COMPILER>.tar` renamed to
`letkf.tar`, `NAMELISTDIR` to `namelists/`, and use `bin/run_LETKF.sh`.
Set `NTEMP` and `ALPHA_S` in `letkf.conf`; if they are unset, the script
writes `ntemp = 1`.
