#!/bin/bash
#PJM -g ra000007
#PJM -L "rscgrp=small"
#PJM -L "node=10"
#PJM -L "elapse=01:00:00"
#PJM --mpi "max-proc-per-node=1"
#PJM -x PJM_LLIO_GFSCACHE=/vol0003:/vol0004
#PJM --llio sharedtmp-size=20Gi
#PJM -j
#PJM -s
#=======================================================================
# ntemp = 1 regression test on Fugaku compute nodes (FUJITSU_FUGAKU build).
#
# Two comparisons, each run by tests/run_regression.sh:
#   A. reference build vs main build      : the change in the source code.
#   B. director's letkf.exe vs reference  : that our build is equivalent to
#                                           the one used in the experiments
#                                           (only needed once).
# Each comparison runs its first executable twice (determinism check).
#
# Before submitting, on the login node:
#   cd ~/WRF_TEnKF
#   ./make.sh FUJITSU_FUGAKU                                  # main
#   git worktree add ../WRF_TEnKF_reference reference
#   (cd ../WRF_TEnKF_reference && ./make.sh FUJITSU_FUGAKU)   # reference
#   tests/make_case.sh $EXP 20240319180000 10 ~/tenkf_tests/case_20240319180000_nbv10
#   pjsub tests/jobs/regression_fugaku.sh
#
# One MPI rank per node, as in the experiments (LETKFSKIP = ICORE = 48),
# with OpenMP threads for speed. NBV members need NBV nodes to read the
# ensemble in one pass; change "node" above together with NBV.
#
# Work files go to the LLIO shared temporary area and vanish with the job;
# logs, NOUT files and obs.dat are copied to RESULTS.
#=======================================================================
set -e
ROOT=${ROOT:-$HOME/WRF_TEnKF}
REF_ROOT=${REF_ROOT:-$HOME/WRF_TEnKF_reference}
CASE_DIR=${CASE_DIR:-$HOME/tenkf_tests/case_20240319180000_nbv10}
DIRECTOR_EXE=${DIRECTOR_EXE:-/data/ra000007/jruiz/data_assimilation_exps/PREVENIR_HIRESOLUTION_DA_FUGAKU_2024031900_5MIN_WINIT/LETKF/code/letkf.exe}
RUN_B=${RUN_B:-1}                       # 0 skips comparison B
NPROC=${NPROC:-10}                      # MPI ranks = nodes
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-48}
export PARALLEL=$OMP_NUM_THREADS
export OMP_STACKSIZE=2G

#Runtime environment of wrfutil (machine.conf, MACHINE=FUGAKU)
export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib/:/vol0004/ra000007/data/jruiz/comp_libs_fujitsu/netcdf/lib:$LD_LIBRARY_PATH
export MPIEXE="mpiexec -n $NPROC"

STAMP=$(date +%Y%m%d_%H%M%S)
RESULTS=$HOME/tenkf_tests/results/regression_${STAMP}_${PJM_JOBID}
WORK=${PJM_SHAREDTMP:?Run this script as a PJM job}/regression
mkdir -p $RESULTS $WORK

{
echo "=== Environment"
echo "job $PJM_JOBID, nodes $PJM_NODE, ranks $NPROC, threads $OMP_NUM_THREADS, host $(hostname)"
echo "main      : $(git -C $ROOT describe --always --dirty --tags) $(md5sum < $ROOT/build/FUJITSU_FUGAKU/letkf.exe)"
echo "reference : $(git -C $REF_ROOT describe --always --dirty --tags) $(md5sum < $REF_ROOT/build/FUJITSU_FUGAKU/letkf.exe)"
echo "director  : $DIRECTOR_EXE $(md5sum < $DIRECTOR_EXE)"
echo "case      : $CASE_DIR"
} | tee $RESULTS/summary.txt

save_logs () {   # $1 = comparison name
  local W=$WORK/$1
  mkdir -p $RESULTS/$1
  cp $W/*.log $W/letkf.namelist $RESULTS/$1/ 2>/dev/null || true
  for r in ref1 ref2 new; do
    [ -d $W/$r ] || continue
    mkdir -p $RESULTS/$1/$r
    cp $W/$r/NOUT-00000 $W/$r/obs.dat $RESULTS/$1/$r/ 2>/dev/null || true
  done
}

status=0
echo "=== A. reference vs main (ntemp = 1)" | tee -a $RESULTS/summary.txt
set +e
$ROOT/tests/run_regression.sh $CASE_DIR $REF_ROOT/build/FUJITSU_FUGAKU/letkf.exe \
     $ROOT/build/FUJITSU_FUGAKU/letkf.exe $WORK/A 2>&1 | tee -a $RESULTS/summary.txt
rcA=${PIPESTATUS[0]}
set -e
save_logs A
[ $rcA -ne 0 ] && status=$rcA

if [ $RUN_B -eq 1 ]; then
  rm -rf $WORK/A                                    # free the shared temporary area
  echo "=== B. director's letkf.exe vs reference build" | tee -a $RESULTS/summary.txt
  set +e
  $ROOT/tests/run_regression.sh $CASE_DIR $DIRECTOR_EXE \
       $REF_ROOT/build/FUJITSU_FUGAKU/letkf.exe $WORK/B 2>&1 | tee -a $RESULTS/summary.txt
  rcB=${PIPESTATUS[0]}
  set -e
  save_logs B
  [ $rcB -ne 0 ] && [ $status -eq 0 ] && status=$rcB
fi

echo "=== Results in $RESULTS (exit $status)" | tee -a $RESULTS/summary.txt
exit $status
