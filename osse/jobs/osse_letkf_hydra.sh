#!/bin/bash
#PBS -N osse_letkf
#PBS -q larga
#PBS -l nodes=12:ppn=48
#PBS -l walltime=06:00:00
#PBS -j oe
#=======================================================================
# OSSE analyses on Hydra compute nodes: osse/osse_letkf.sh once for each
# number of tempering steps in NTEMP_LIST, in one job. The Hydra
# counterpart of osse_letkf_fugaku.sh.
#
# One MPI process per node with 48 OpenMP threads (osse_letkf.sh builds
# the machine file from the job's nodes). The ensemble (59 members, prior
# and analysis, about 310 GB in double precision) is spread over the
# nodes: 12 nodes of 128 GB leave room for the observations.
#
# Submit from the repository, in the scratch area:
#   cd ~/salidas/WRF_TEnKF
#   qsub osse/jobs/osse_letkf_hydra.sh                                 # defaults below
#   qsub -v NTEMP_LIST=1:2:4,SKIP_IJ=2,SKIP_K=2 osse/jobs/osse_letkf_hydra.sh
#
# Variables (qsub -v NAME=value,...; defaults in brackets):
#   GUES_DIR   [$BASE/osse/da4km/HIST/GUES/20240319180000]
#   NMEM [60]  TRUTH [60]
#   NTEMP_LIST ["1:2"], separated by ':'
#   ALPHA_S [1.0d0]  TEMP_RELAX_MODE [1]
#   SKIP_IJ, SKIP_K [1, 1]  thinning of the synthetic observations (the
#              Fugaku OSSE used 2, 2: 669k observations at 4 km)
#   OMP_NUM_THREADS [48]
#   OUT_ROOT   [$BASE/osse/osse18]: results in $OUT_ROOT/anal_nt<N> (an
#              existing directory is skipped, so a resubmission only runs
#              what is missing)
# Each analysis copies the NMEM-1 guesses to a work folder that it deletes
# at the end.
#=======================================================================
BASE=${BASE:-$HOME/salidas}
ROOT=${ROOT:-$BASE/WRF_TEnKF}
GUES_DIR=${GUES_DIR:-$BASE/osse/da4km/HIST/GUES/20240319180000}
NMEM=${NMEM:-60}
TRUTH=${TRUTH:-60}
NTEMP_LIST=${NTEMP_LIST:-1:2}
NTEMP_LIST=${NTEMP_LIST//:/ }
ALPHA_S=${ALPHA_S:-1.0d0}
TEMP_RELAX_MODE=${TEMP_RELAX_MODE:-1}
OUT_ROOT=${OUT_ROOT:-$BASE/osse/osse18}
export SKIP_IJ=${SKIP_IJ:-1} SKIP_K=${SKIP_K:-1}
[ -n "$PBS_NODEFILE" ] || { echo "Run this script as a PBS job"; exit 1; }

echo "== job ${PBS_JOBID}, nodes $(sort -u $PBS_NODEFILE | tr '\n' ' ')"
echo "   guesses $GUES_DIR, truth $TRUTH, ntemp in: $NTEMP_LIST"
mkdir -p $OUT_ROOT
RC=0
for NT in $NTEMP_LIST; do
  OUT=$OUT_ROOT/anal_nt$NT
  if [ -e $OUT ]; then echo "== $OUT exists, skipped"; continue; fi
  if $ROOT/osse/osse_letkf.sh $GUES_DIR $NMEM $TRUTH $OUT $NT $ALPHA_S $TEMP_RELAX_MODE; then
    echo "== ntemp=$NT done"
  else
    echo "== ntemp=$NT FAILED"; RC=1
    mv $OUT $OUT.failed.$(date +%s) 2> /dev/null   #so that a resubmission runs it again
  fi
done
exit $RC
