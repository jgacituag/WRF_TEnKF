#!/bin/bash
#PJM -g ra000007
#PJM -L "rscgrp=small"
#PJM -L "node=60"
#PJM -L "elapse=02:00:00"
#PJM --mpi "max-proc-per-node=1"
#PJM -x PJM_LLIO_GFSCACHE=/vol0004
#PJM -j
#PJM -s
#=======================================================================
# OSSE analyses on Fugaku compute nodes: osse/osse_letkf.sh once for each
# number of tempering steps in NTEMP_LIST, in one job (one queue wait).
# One MPI rank per node with 48 OpenMP threads, as in the experiments;
# with at least NMEM-1 nodes the LETKF reads the ensemble in one pass.
#
# Submit from the repository, in the data area:
#   cd /data/ra000007/$USER/WRF_TEnKF
#   pjsub osse/jobs/osse_letkf_fugaku.sh                       # defaults below
#   pjsub -x NTEMP_LIST=1:2:4,ALPHA_S=1.0d0 osse/jobs/osse_letkf_fugaku.sh
#
# Variables (pjsub -x NAME=value,...; defaults in brackets):
#   GUES_DIR   [$BASE/osse/da4km/HIST/GUES/20240319180000]
#   NMEM [60]  TRUTH [60]
#   NTEMP_LIST ["1:2"], separated by ':' (pjsub -x splits on commas)
#   ALPHA_S [1.0d0]  TEMP_RELAX_MODE [1]
#   OUT_ROOT   [$BASE/osse/osse18]: results in $OUT_ROOT/anal_nt<N>
#              (an existing directory is skipped, so a resubmission only
#              runs what is missing)
# Each analysis copies the NMEM-1 guesses (0.7 GB each) to a work folder
# that it deletes at the end.
#=======================================================================
BASE=${BASE:-/data/ra000007/$USER}
ROOT=${ROOT:-$BASE/WRF_TEnKF}
GUES_DIR=${GUES_DIR:-$BASE/osse/da4km/HIST/GUES/20240319180000}
NMEM=${NMEM:-60}
TRUTH=${TRUTH:-60}
NTEMP_LIST=${NTEMP_LIST:-1:2}
NTEMP_LIST=${NTEMP_LIST//:/ }
ALPHA_S=${ALPHA_S:-1.0d0}
TEMP_RELAX_MODE=${TEMP_RELAX_MODE:-1}
OUT_ROOT=${OUT_ROOT:-$BASE/osse/osse18}
export NPROC=${NPROC:-${PJM_NODE:-$((NMEM-1))}}

echo "== job ${PJM_JOBID:-local}, $NPROC ranks, guesses $GUES_DIR, truth $TRUTH, ntemp in: $NTEMP_LIST"
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
