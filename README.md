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

## Use with wrfutil

`wrfutil_hydra` stays unchanged. In the experiment configuration, point
`LETKFDIR` to a directory holding `dist/letkf_<COMPILER>.tar` renamed to
`letkf.tar`, `NAMELISTDIR` to `namelists/`, and use `bin/run_LETKF.sh`.
Set `NTEMP` and `ALPHA_S` in `letkf.conf`; if they are unset, the script
writes `ntemp = 1`.
