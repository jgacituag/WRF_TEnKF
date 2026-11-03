#!/bin/bash
#=======================================================================
# Unit tests for the tempering weights and the &tempering namelist block.
# Needs a previous ./make.sh <COMPILATION_NAME> (default GNU_LOCAL), and
# runs serially, without MPI or WRF files.
#=======================================================================
set -e
COMP=${1:-GNU_LOCAL}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
BUILD=$ROOT/build/$COMP
[ -d $BUILD ] || { echo "Run ./make.sh $COMP first"; exit 1; }
case $COMP in
  GNU_LOCAL) F90="mpif90 -fopenmp -ffree-line-length-none"; INC="$(nf-config --fflags)"; LIB="$(nf-config --flibs)" ;;
  *) echo "Set the compiler for $COMP in this script"; exit 1 ;;
esac

WORK=$(mktemp -d)
trap "rm -rf $WORK" EXIT
cd $BUILD
OBJS="SFMT.o common.o module_map_utils.o common_wrf.o common_namelist.o"
$F90 $INC -o $WORK/test_weights.exe $ROOT/tests/test_weights.f90 $OBJS $LIB
$F90 $INC -o $WORK/test_namelist_read.exe $ROOT/tests/test_namelist_read.f90 $OBJS $LIB
cd $WORK
nfail=0

echo "== Weights against thesis Table 2.1"
./test_weights.exe || nfail=$((nfail+1))

echo "== Namelist without &tempering: expect ntemp=1"
printf '&general\n nbv=60\n/\n' > letkf.namelist
out=$(./test_namelist_read.exe | grep ntemp=)
echo "$out"; [[ "$out" == *"ntemp=  1 temp_alpha=  1.000000"* ]] || { echo FAIL; nfail=$((nfail+1)); }

echo "== Namelist with &tempering ntemp=3 alpha_s=1"
printf '&general\n nbv=60\n/\n&tempering\n ntemp=3\n alpha_s=1.0d0\n/\n' > letkf.namelist
out=$(./test_namelist_read.exe | grep ntemp=)
echo "$out"; [[ "$out" == *"ntemp=  3 temp_alpha=  0.043896  0.324352  0.631752"* ]] || { echo FAIL; nfail=$((nfail+1)); }

echo "== Namelist with a misspelled &tempering entry: expect a stop"
printf '&tempering\n ntmep=3\n/\n' > letkf.namelist
if ./test_namelist_read.exe > log 2>&1 ; then echo "FAIL: did not stop"; nfail=$((nfail+1)); else grep ERROR log; fi

echo "== alpha_list not summing to one: expect a stop"
printf '&tempering\n ntemp=2\n alpha_list=0.5d0,0.6d0\n/\n' > letkf.namelist
if ./test_namelist_read.exe > log 2>&1 ; then echo "FAIL: did not stop"; nfail=$((nfail+1)); else grep ERROR log; fi

[ $nfail -eq 0 ] && echo "ALL TESTS PASSED" || { echo "$nfail TEST GROUPS FAILED"; exit 1; }
