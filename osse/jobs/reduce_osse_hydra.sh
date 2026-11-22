#!/bin/bash
#PBS -N reduce_osse
#PBS -q larga
#PBS -l nodes=1:ppn=48
#PBS -l walltime=03:00:00
#PBS -j oe
#=======================================================================
# Reduce one OSSE run on one Hydra compute node: verif/reduce_osse.py for
# EXP (noda or ntN) at TIMES, reading the run from the scratch area and
# writing to OUT. Reduce each forecast run before it leaves the scratch
# area. The compute nodes do not reach the datosmunin disks: copy OUT there
# from the login node afterwards, for the notebooks.
#
# Submit from the repository, in the scratch area, one job per run (it
# takes a whole node: 48 cores, 128 GB):
#   cd ~/salidas/WRF_TEnKF
#   qsub -v EXP=nt1,TIMES=1800 osse/jobs/reduce_osse_hydra.sh   # OSSE analysis only
#   qsub -v EXP=noda osse/jobs/reduce_osse_hydra.sh              # prior at 18:00, forecast to 19:00
#
# Variables (qsub -v NAME=value,...; defaults in brackets):
#   EXP     run to reduce (required): noda or ntN
#   TIMES   valid times HHMM separated by ':' [18:00 to 19:00 every 5 min]
#   OUT     output folder [$BASE/osse/osse18/verif]
#   NPROC   reading processes [16]
#   PY      Python with numpy, scipy and netCDF4 [$BASE/miniconda3/bin/python3]
#=======================================================================
BASE=${BASE:-$HOME/salidas}
ROOT=${ROOT:-$BASE/WRF_TEnKF}
PY=${PY:-$BASE/miniconda3/bin/python3}
OUT=${OUT:-$BASE/osse/osse18/verif}
NPROC=${NPROC:-16}
[ -n "$EXP" ] || { echo "EXP not set: qsub -v EXP=ntN[,TIMES=1800] osse/jobs/reduce_osse_hydra.sh"; exit 1; }
[ -x "$PY" ] || { echo "no Python at $PY"; exit 1; }

ARGS="--osse $BASE/osse --exp $EXP --out $OUT --nproc $NPROC"
[ -n "$TIMES" ] && ARGS="$ARGS --times ${TIMES//:/,}"
export HDF5_USE_FILE_LOCKING=FALSE
echo "== job ${PBS_JOBID} on $(hostname): $PY verif/reduce_osse.py $ARGS"
cd $ROOT && $PY -u verif/reduce_osse.py $ARGS
