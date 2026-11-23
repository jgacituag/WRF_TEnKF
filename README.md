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

How each step works (`src/letkf/letkf.f90`):

- Step 1 is the original LETKF with R/α_1: `set_letkf_obs` reads the
  observations, applies the quality control and H(x).
- Steps 2..Nt keep the observations, their quality control and their
  order from step 1. `set_letkf_obs_temp` gathers the previous analysis
  member by member in memory (`gather_ens_mpi_alltoall`, the
  communication of `write_ens_mpi_alltoall` without the files) and
  recomputes H(x) and the departures. An observation the operator cannot
  evaluate at a step (it left the vertical domain, say) has no effect in
  that step.
- Moisture and hydrometeors are clipped at zero after every step, and MU
  is updated against that step's prior, as in the original code.
- Fixed multiplicative inflation ρ is split among the steps: step i uses
  ρ^α_i, so the steps together inflate by ρ (the same for
  `min_infl_mul`). Adaptive inflation is not supported with Nt > 1.
- RTPP/RTPS: `temp_relax_mode = 1` (default) relaxes once, after the last
  step, against the original prior, in state space; `temp_relax_mode = 2`
  relaxes at every step against that step's prior (the original code in
  weight space).
- Nt > 1 requires one time slot (`nslots = 1`) and no parameter
  estimation.

## Namelist

```
&tempering
  ntemp   = 1        ! 1 = original LETKF
  alpha_s = 1.0d0    ! delay parameter of the exponential weights
  alpha_list = ...   ! optional explicit weights (override alpha_s)
  temp_relax_mode = 1  ! RTPP/RTPS: 1 = once after the last step, 2 = every step
/
```

A missing `&tempering` block leaves `ntemp = 1`.

## Layout

```
src/common/      generic LETKF modules (Miyoshi)
src/wrf_common/  WRF interface, MPI I/O, observation operators, namelist
src/letkf/       main program, analysis (letkf_tools), observations (letkf_obs)
src/osse/        synth_obs.f90, synthetic reflectivity from a truth member
namelists/       letkf.namelist template for bin/run_LETKF.sh; OSSE settings
bin/             run_LETKF.sh for wrfutil, with the tempering placeholders
osse/            experiment builders (DA cycle, forecasts) and the OSSE analysis
osse/jobs/       batch jobs: OSSE analyses (PJM, PBS), reduction (PBS)
osse/conf/       wrfutil template of the 4 km, 5-min DA cycle on Hydra
verif/           OSSE verification in Python
tests/           ntemp = 1 regression test and weights test
tests/jobs/      regression jobs for Fugaku (PJM) and Hydra (PBS)
tests/local/     workstation tests on small synthetic files
make.sh          build script; writes dist/letkf_<COMPILER>.tar
```

## Running an experiment

The steps below run an OSSE from a DA cycle: build, test, cycle, OSSE
analyses, forecasts and verification. Commands run from the repository
root unless a step says otherwise.

Machine requirements:

- **Fugaku:** keep the repository, builds, test cases and results in the
  group data area, not in the home area (20 GiB per user): jobs must be
  submitted from the data area and read their files from it. Directories
  that will receive files larger than 1 GB from compute nodes (wrfout)
  need striping, `lfs setstripe -c 8 <dir>`.
- **Hydra (CIMA):** keep everything in the scratch area, which the
  compute nodes reach. `letkf.exe` links Intel MKL from the login
  environment (oneAPI 2026.0); the Hydra jobs and the cycle template
  (`machine.conf`) add its folder to `LD_LIBRARY_PATH`. Disks outside the scratch area that only the login
  node reaches (archives) are read and written from the login node only.

### 1. Build

```
./make.sh INTEL_HYDRA                  # or FUJITSU_FUGAKU, INTEL_FUGAKU
git worktree add ../WRF_TEnKF_reference reference
(cd ../WRF_TEnKF_reference && ./make.sh INTEL_HYDRA)   # the original code, its own make.sh
```

Each build writes `dist/letkf_<COMPILER>.tar` (`letkf.exe`,
`synth_obs.exe` and the other executables).

### 2. Test

```
tests/test_weights.sh INTEL_HYDRA
tests/make_case.sh EXPDIR DATE 10 CASE_DIR
qsub -v CASE_DIR=CASE_DIR tests/jobs/regression_hydra.sh   # Hydra
pjsub tests/jobs/regression_fugaku.sh                      # Fugaku (settings in its header)
```

`make_case.sh` needs a wrfutil experiment that has run at least one
analysis at `DATE`. The regression test must report bit-for-bit
identical outputs before any experiment (see Tests).

### 3. DA cycle

`wrfutil_hydra` runs the cycle and stays unchanged. In the experiment
configuration, point `LETKFDIR` to a directory holding
`dist/letkf_<COMPILER>.tar` renamed to `letkf.tar`, `NAMELISTDIR` to
`namelists/`, and use `bin/run_LETKF.sh`. Set `NTEMP`, `ALPHA_S` and
`TEMP_RELAX_MODE` in `letkf.conf`; if they are unset, the script writes
`ntemp = 1`.

- Hydra: `osse/make_exp_hydra.sh [--dry-run] NEW_EXP ["DA_INI_DATE"
  "DA_END_DATE"]` builds the 4 km, 5-min cycle from `osse/conf/` with the
  `reference` build (`--letkf` chooses another tar; `--smoke` makes a
  4-member test).
- Fugaku: `osse/make_da_exp.sh [--dry-run] SRC_EXP NEW_EXP "DA_INI_DATE"
  "DA_END_DATE"` repeats an existing experiment over a new period.

Then `cd NEW_EXP/bin && ./queue_sub.sh main_ASIM.sh`. After a time limit,
set `INI_STEP` in `NEW_EXP/conf/exp.conf` to the first unfinished step and
submit again.

### 4. OSSE analyses

```
qsub -v NTEMP_LIST=1:2:3:4:5,SKIP_IJ=2,SKIP_K=2 osse/jobs/osse_letkf_hydra.sh   # Hydra
pjsub -x NTEMP_LIST=1:2:4 osse/jobs/osse_letkf_fugaku.sh                         # Fugaku
```

One analysis per value of `NTEMP_LIST`, from the guesses in `GUES_DIR`
with member `TRUTH` as the truth; each writes `OUT_ROOT/anal_nt<N>`. The
variables and their defaults are listed in each job's header.

### 5. Forecasts

```
osse/link_gues_as_anal.sh GUES_DIR 60 ANAL_NODA
osse/make_fcst_exp.sh --dry-run --gefs DA_EXP FC_NODA "YYYY-MM-DD hh:mm:ss" 3600 001 060 ANAL_NODA
osse/make_fcst_exp.sh --gefs DA_EXP FC_NODA "YYYY-MM-DD hh:mm:ss" 3600 001 060 ANAL_NODA
osse/make_fcst_exp.sh FC_NODA FC_NT2 "YYYY-MM-DD hh:mm:ss" 3600 001 060 OUT_ROOT/anal_nt2
cd FC_NODA/bin && ./queue_sub.sh main_FCST.sh
```

The forecast without DA (`--gefs`) runs WPS and the boundary
perturbations as the DA cycle did; it also forecasts the truth member.
The forecasts from the OSSE analyses link its `met_em`, so keep
`FC_NODA/HIST/WPS` until all of them have run. Outputs every 5 min go to
`HIST/DAFCST/<date>/<member>/`; after each run, its `clean_work.sh`
removes the work folders.

### 6. Verify

```
qsub -v EXP=noda osse/jobs/reduce_osse_hydra.sh
qsub -v EXP=nt2 osse/jobs/reduce_osse_hydra.sh
```

One job per run, as soon as the run ends (`TIMES` limits the valid
times). Every reduction reads the truth: the truth member's guess at the
analysis time and its forecast in `FC_NODA` after it, so keep both until
every run is reduced. Copy the reduced files (`OUT`) to where the
notebooks read them (`VERIF_DIR`) and run the notebooks from there. Once,
before the first reduction, `verif/check_operator.py` checks the Python
reflectivity operator against `synth_obs` (see OSSE verification).

## Tests

- `tests/test_weights.sh [COMPILER]` runs serially, without WRF data. It
  checks the weights against thesis Table 2.1 and the reading of the
  `&tempering` block (absent, present, misspelled, weights not summing
  to one).
- `tests/run_regression.sh CASE_DIR REF_EXE NEW_EXE WORK_DIR` runs inside
  one job. It runs the original executable twice and this one once with
  `ntemp = 1`, on fresh copies of the same inputs, and compares every
  file `letkf.exe` writes with `cmp`. Exit code 0 means bit-for-bit
  identical, 1 means different, 2 means the original code is not
  reproducible in that layout. Run it after every change to `src/` and
  before every experiment.
- `tests/make_case.sh EXPDIR DATE NBV CASE_DIR` builds a test case from a
  saved wrfutil experiment, with links only. wrfutil keeps no prior
  members, so the case uses the saved analyses as "prior": meaningless
  physically, identical for every executable.
- `tests/jobs/regression_fugaku.sh` (PJM) compares the `reference` build
  against `main` (the source change) and, once, the `letkf.exe` of the
  experiments against the `reference` build (that our build is
  equivalent). `tests/jobs/regression_hydra.sh` (PBS) makes the first
  comparison only: on Hydra the DA cycle itself runs the `reference`
  build. Instructions are in their headers.

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
- `osse/make_exp_hydra.sh` builds the 4 km, 5-min DA cycle on Hydra from
  the template `osse/conf/PREVENIR_4KM_5MIN_HYDRA_2024031900`: the DA
  settings of the director's Fugaku run (`..._5MIN_WINIT`) with the Hydra
  machine and multi-physics model settings of wrfutil's
  `PREVENIR_LOWRESOLUTION_DA_HYDRA`. The wrfutil clone is not modified.
- `osse/make_da_exp.sh` builds a wrfutil DA-cycle experiment that repeats
  an existing one over a new period (same scripts, namelists, executables,
  GEFS boundaries and radar data; only the dates, `BASEDIR` and
  `EXPTYPE=DACYCLE` change). `--dry-run` checks the executables, the GEFS
  files of every step and the radar files; `--smoke` makes a 4-member,
  8-node version for a check in the interactive group.
- `osse/make_fcst_exp.sh` builds a wrfutil DAFCST experiment that
  forecasts from a given set of analyses, reusing the scripts, namelists,
  WRF executables and met_em of an existing experiment (copied or linked,
  never modified), on Fugaku and on Hydra (defaults follow the source
  experiment's `machine.conf`). `--dry-run` checks every input first.
- `osse/osse_letkf.sh` runs one OSSE analysis inside a job: copies the
  guesses of a DA run (never writes into it), synthetic reflectivity from
  the truth member, LETKF with the requested tempering; writes the
  analyses under the original member numbers plus a copy of the truth
  guess. `osse/jobs/osse_letkf_fugaku.sh` (PJM) and
  `osse/jobs/osse_letkf_hydra.sh` (PBS) run it for a list of `ntemp`
  values. `osse/link_gues_as_anal.sh` prepares the forecast without DA.
- `tests/local/test_osse_tools.sh` runs the chain on small synthetic
  wrfout files on a workstation (needs numpy and netCDF4; nothing of it
  runs on the cluster).
- `tests/local/test_tempering.sh [BASE]` checks the tempering loop on the
  same synthetic files: (A) `ntemp = 1` bit for bit against revision BASE,
  with and without inflation, RTPP and RTPS; (B) with temperature
  observations, global localization and no inflation, where tempering is
  exact algebra, `ntemp = 2` and `4` reproduce the single-step analysis
  (to single precision); (C) reflectivity with `ntemp = 3` in every RTPP
  mode.

## OSSE verification (`verif/`)

Python 3 with numpy, scipy, pandas, matplotlib and netCDF4.

- `osse_io.py` reads wrfout files and applies the reflectivity operator of the LETKF
  (Tong & Xue 2006, method 2 of `calc_ref_vr`) in NumPy.
- `check_operator.py TRUTH_WRFOUT truth_rad0101.dat` compares it with the Fortran
  values `synth_obs` wrote for the truth at every observation (must agree within
  0.01 dBZ).
- `reduce_osse.py` reduces one run (`--exp noda` or `ntN`) at the requested valid
  times, reading only that run's members and the truth, so that each forecast run can be
  reduced as soon as it ends. It writes pooled sums per variable over domains fixed by
  the truth (whole box, observed cells, storm columns, truth echo), in total and per
  level, and compact per-cell fields: reflectivity mean, variance and CRPS, the mean
  hydrometeor mass and the column maximum of every member; plus the truth at each time
  and the masks. `osse/jobs/reduce_osse_hydra.sh` runs it on one node of Hydra.
- `vcommon.py` holds what the verification notebooks share: loading the reduced files;
  RMSE, bias, CRPS, spread, spread/RMSE and skill scores; what compares runs at one
  time (the echo domain shared by all runs, false and missed echo and the hydrometeor
  mass outside the truth's echo, the fraction of cells that improve, rank histograms and
  probabilities of the column maximum); FSS, ROC AUC and frequency bias of the column
  maximum (`spatial_scores`); and bilingual figures. The notebooks are kept outside the
  repository, next to the reduced files (`VERIF_DIR`).
