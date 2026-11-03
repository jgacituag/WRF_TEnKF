#!/bin/bash
#=======================================================================
# Regression test: WRF_TEnKF with ntemp = 1 must reproduce the original
# WRF-LETKF bit for bit.
#
# Runs three analyses on the same inputs, in the same job (same nodes,
# MPI ranks and OpenMP threads):
#   ref1, ref2 : original letkf.exe (wrfutil_hydra@586cc96), twice, to
#                confirm that the original code is deterministic here;
#   new        : WRF_TEnKF letkf.exe with ntemp = 1.
# Then it compares the files letkf.exe writes: gs<NBSLOT><member>,
# guesemean, analemean and obs.dat.
#
# letkf.exe overwrites gs<NBSLOT>* in place, so each run gets its own
# fresh copy of those files. Everything else is linked.
#
# Usage (inside a job script):
#   export MPIEXE="mpiexec"            # or "mpiexec -n 60", "srun", ...
#   export OMP_NUM_THREADS=...
#   tests/run_regression.sh CASE_DIR REF_EXE NEW_EXE WORK_DIR
#
# CASE_DIR holds the inputs with the names run_LETKF.sh gives them in
# $LETKFDIR/00/ before calling letkf.exe: letkf.namelist (placeholders
# already replaced), gs<slot><member>, obs<slot>.dat, rad<slot><radar>.dat,
# guesemean and analemean. They may be links, but they must point to
# PRIOR files (forecast wrfouts): in a normal cycle letkf.exe overwrites
# the wrfouts that gs* link to with the analysis. This script never
# writes into CASE_DIR.
#=======================================================================
set -e
if [ $# -ne 4 ]; then
  sed -n '3,/^#====/p' "$0"; exit 1
fi
CASE_DIR=$(cd "$1" && pwd)
REF_EXE=$(readlink -f "$2")
NEW_EXE=$(readlink -f "$3")
WORK_DIR=$4
MPIEXE=${MPIEXE:?Set MPIEXE, e.g. export MPIEXE=mpiexec}
HERE=$(cd "$(dirname "$0")" && pwd)

[ -f $CASE_DIR/letkf.namelist ] || { echo "No letkf.namelist in $CASE_DIR"; exit 1; }
NBSLOT=$(grep -i "^ *nbslot *=" $CASE_DIR/letkf.namelist | sed 's/.*= *\([0-9]*\).*/\1/')
[ -n "$NBSLOT" ] || { echo "Could not read nbslot from letkf.namelist"; exit 1; }
WRITTEN="$(cd $CASE_DIR && ls gs$(printf %02d $NBSLOT)?????) guesemean analemean"

#The namelist must ask for ntemp = 1. Without a &tempering block, add one
#so that the new code also exercises the namelist reading path.
NML=$WORK_DIR/letkf.namelist
mkdir -p $WORK_DIR
cp $CASE_DIR/letkf.namelist $NML
if grep -qi "&tempering" $NML ; then
  NT=$(sed -n '/&tempering/I,/\//p' $NML | grep -i "ntemp *=" | sed 's/.*= *\([0-9]*\).*/\1/')
  [ "$NT" == "1" ] || { echo "letkf.namelist has ntemp=$NT; the regression test needs ntemp=1"; exit 1; }
else
  printf '&tempering\n  ntemp = 1\n  alpha_s = 1.0d0\n/\n' >> $NML
fi

run_case () {   # $1 = run name, $2 = executable
  local RUN=$WORK_DIR/$1
  rm -rf $RUN && mkdir -p $RUN
  for f in $(cd $CASE_DIR && ls); do
    [ "$f" == "letkf.namelist" ] && continue
    ln -s $(readlink -f $CASE_DIR/$f) $RUN/$f
  done
  for f in $WRITTEN; do                # fresh, real copies of the files letkf.exe writes
    rm -f $RUN/$f && cp -L $CASE_DIR/$f $RUN/$f
  done
  cp $NML $RUN/letkf.namelist
  ln -s $2 $RUN/letkf.exe
  echo "== Running $1 ($2)"
  ( cd $RUN && export FORT90L="" && time $MPIEXE ./letkf.exe ) > $WORK_DIR/$1.log 2>&1
  grep -q "LETKF ends SUCCESFULLY" $RUN/NOUT-00000 || { echo "$1 did not finish, see $RUN/NOUT-00000"; exit 1; }
}

compare () {    # $1, $2 = run names. Returns the number of differing files.
  local nd=0
  for f in $WRITTEN obs.dat; do
    if ! cmp -s $WORK_DIR/$1/$f $WORK_DIR/$2/$f ; then
      echo "  DIFFERENT: $f"; nd=$((nd+1))
    fi
  done
  if [ $nd -gt 0 ] && command -v python3 > /dev/null ; then
    python3 $HERE/compare_outputs.py $WORK_DIR/$1 $WORK_DIR/$2 $WRITTEN || true
  fi
  return $nd
}

run_case ref1 $REF_EXE
run_case ref2 $REF_EXE
run_case new  $NEW_EXE
grep "ntemp\|temp_alpha" $WORK_DIR/new/NOUT-00000 || \
  echo "WARNING: the new NOUT does not report ntemp; is NEW_EXE the WRF_TEnKF build?"

echo "== ref1 vs ref2 (determinism of the original code)"
if compare ref1 ref2 ; then
  echo "  identical"
else
  echo "INCONCLUSIVE: the original code is not reproducible in this setup."
  echo "Fix the layout (ranks, threads, nodes) until ref1 == ref2."
  exit 2
fi

echo "== ref1 vs new (ntemp = 1)"
if compare ref1 new ; then
  echo "PASSED: ntemp = 1 reproduces the original LETKF bit for bit ($(echo $WRITTEN | wc -w) files + obs.dat)"
else
  echo "FAILED: ntemp = 1 differs from the original LETKF"
  exit 1
fi
