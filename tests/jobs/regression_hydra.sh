#!/bin/bash
#PBS -N tenkf_regression
#PBS -q larga
#PBS -l nodes=10:ppn=48
#PBS -l walltime=02:00:00
#PBS -j oe
#=======================================================================
# ntemp = 1 regression test on Hydra compute nodes (INTEL_HYDRA build):
# the reference build (original code) against the main build, with
# tests/run_regression.sh, which runs the reference twice (determinism
# check) and compares every file letkf.exe writes. The Hydra counterpart of
# regression_fugaku.sh, comparison A only: on Hydra the DA cycle itself
# runs the reference build, so there is no third executable to compare.
#
# Everything lives in the scratch area, BASE = ~/salidas. Before
# submitting, on the login node:
#   cd $BASE/WRF_TEnKF && ./make.sh INTEL_HYDRA                      # main
#   git worktree add ../WRF_TEnKF_reference reference
#   (cd ../WRF_TEnKF_reference && ./make.sh INTEL_HYDRA)            # its own make.sh
#   tests/make_case.sh $EXP <DATE> 10 $BASE/tenkf_tests/case_<DATE>_nbv10
#   qsub -v CASE_DIR=$BASE/tenkf_tests/case_<DATE>_nbv10 tests/jobs/regression_hydra.sh
# with EXP a wrfutil experiment that has run at least one analysis at
# DATE (HIST/ANAL/DATE and LETKF/00/letkf.namelist), e.g. the DA cycle.
#
# One MPI process per node with 48 OpenMP threads, as the DA cycle
# (LETKFSKIP = ICORE = 48, LETKFOMP = 1). NBV members need NBV nodes to
# read the ensemble in one pass; change "nodes" above together with NBV.
#
# Work files go to $BASE/tenkf_tests/work_<job> and are deleted at the
# end; logs, NOUT files and obs.dat are copied to RESULTS.
#=======================================================================
set -e
BASE=${BASE:-$HOME/salidas}
ROOT=${ROOT:-$BASE/WRF_TEnKF}
REF_ROOT=${REF_ROOT:-$BASE/WRF_TEnKF_reference}
CASE_DIR=${CASE_DIR:?Set CASE_DIR (qsub -v CASE_DIR=...)}
BUILD=INTEL_HYDRA

#Runtime environment of wrfutil (machine.conf, MACHINE=HYDRA)
source /opt/load-libs.sh 1
export HDF5_USE_FILE_LOCKING=FALSE
ulimit -s unlimited
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-48}
export OMP_STACKSIZE=2G

JOBID=${PBS_JOBID%%.*}
STAMP=$(date +%Y%m%d_%H%M%S)
RESULTS=$BASE/tenkf_tests/results/regression_${STAMP}_${JOBID}
WORK=$BASE/tenkf_tests/work_${JOBID:-local}
mkdir -p $RESULTS $WORK
MACHINES=$WORK/machines
sort -u ${PBS_NODEFILE:?Run this script as a PBS job} > $MACHINES
NPROC=$(wc -l < $MACHINES)
export MPIEXE="/opt/intel/oneapi/mpi/latest/bin/mpiexec -np $NPROC -machinefile $MACHINES"

{
echo "=== Environment"
echo "job $PBS_JOBID, nodes $(tr '\n' ' ' < $MACHINES), processes $NPROC, threads $OMP_NUM_THREADS, host $(hostname)"
echo "main      : $(git -C $ROOT describe --always --dirty --tags) $(md5sum < $ROOT/build/$BUILD/letkf.exe)"
echo "reference : $(git -C $REF_ROOT describe --always --dirty --tags) $(md5sum < $REF_ROOT/build/$BUILD/letkf.exe)"
echo "case      : $CASE_DIR"
} | tee $RESULTS/summary.txt

echo "=== A. reference vs main (ntemp = 1)" | tee -a $RESULTS/summary.txt
set +e
$ROOT/tests/run_regression.sh $CASE_DIR $REF_ROOT/build/$BUILD/letkf.exe \
     $ROOT/build/$BUILD/letkf.exe $WORK/A 2>&1 | tee -a $RESULTS/summary.txt
status=${PIPESTATUS[0]}
set -e
mkdir -p $RESULTS/A
cp $WORK/A/*.log $WORK/A/letkf.namelist $RESULTS/A/ 2>/dev/null || true
for r in ref1 ref2 new; do
  [ -d $WORK/A/$r ] || continue
  mkdir -p $RESULTS/A/$r
  cp $WORK/A/$r/NOUT-00000 $WORK/A/$r/obs.dat $RESULTS/A/$r/ 2>/dev/null || true
done
rm -rf $WORK
echo "=== Results in $RESULTS (exit $status)" | tee -a $RESULTS/summary.txt
exit $status
