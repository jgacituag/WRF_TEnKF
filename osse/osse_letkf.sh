#!/bin/bash
#=======================================================================
# One OSSE analysis: synthetic reflectivity from a truth member, then the
# LETKF (tempered or not) on the other members. Run it inside a job (e.g.
# an interactive session with one MPI rank per node, like the regression
# test); it reads the guesses of a wrfutil DA experiment and never writes
# into that experiment.
#
# Usage:
#   osse/osse_letkf.sh GUES_DIR NMEM TRUTH OUT_DIR [NTEMP [ALPHA_S [TEMP_RELAX_MODE]]]
#
#   GUES_DIR  directory with gues<NNNNN> files (HIST/GUES/<time> of the DA run)
#   NMEM      members in GUES_DIR (60 in the OSSE, 4 in the smoke test)
#   TRUTH     member used as the truth (060); the other NMEM-1 members are
#             the ensemble
#   OUT_DIR   new directory for the results (must not exist):
#               anal<NNNNN>     analyses, under the original member numbers
#               anal<TRUTH>     copy of the truth guess, so that a forecast
#                               experiment of all NMEM members also runs the
#                               truth forecast (make_fcst_exp.sh --gefs)
#               guesemean, analemean, NOUT-00000, obs.dat, letkf.namelist
#               rad0101.dat, truth_rad0101.dat, SYNO-00000 (synth_obs summary)
#   NTEMP, ALPHA_S, TEMP_RELAX_MODE  tempering (defaults 1, 1.0d0, 1)
#
# Environment (defaults in brackets):
#   NPROC [8]  MPI ranks;  MPIEXE ["mpiexec -n $NPROC"]
#   LETKF_EXE, SYNTH_EXE [build/FUJITSU_FUGAKU of this repository]
#   NAMELIST [namelists/letkf_osse_cap4.namelist]
#   SEED [60] noise seed of the synthetic observations
#   KEEP_WORK [0] keep OUT_DIR/work (copies of the guesses) afterwards
#
# Example, dress rehearsal on the smoke test (8-node interactive session,
# --mpi "max-proc-per-node=1"):
#   osse/osse_letkf.sh $D/osse/smoke4km/HIST/GUES/20240319120500 4 4 \
#                      $D/osse/rehearsal/anal_nt2 2
#=======================================================================
set -e
if [ $# -lt 4 ]; then sed -n '3,/^#====/p' "$0"; exit 1; fi
ROOT=$(cd "$(dirname "$0")/.." && pwd)
GUES=$(readlink -f "$1")
NMEM=$((10#$2))
TRUTH=$((10#$3))
OUT=$(readlink -m "$4")
NTEMP=${5:-1}
ALPHA_S=${6:-1.0d0}
RELAX=${7:-1}
NPROC=${NPROC:-8}
MPIEXE=${MPIEXE:-"mpiexec -n $NPROC"}
LETKF_EXE=${LETKF_EXE:-$ROOT/build/FUJITSU_FUGAKU/letkf.exe}
SYNTH_EXE=${SYNTH_EXE:-$ROOT/build/FUJITSU_FUGAKU/synth_obs.exe}
NAMELIST=${NAMELIST:-$ROOT/namelists/letkf_osse_cap4.namelist}
LETKF_EXE=$(readlink -f $LETKF_EXE) ; SYNTH_EXE=$(readlink -f $SYNTH_EXE) ; NAMELIST=$(readlink -f $NAMELIST)
SEED=${SEED:-60}
NBV=$((NMEM - 1))

#Runtime environment of wrfutil on Fugaku (machine.conf, MACHINE=FUGAKU)
if [ -d /opt/FJSVxtclanga ]; then
  export LD_LIBRARY_PATH=/lib64:/usr/lib64:/opt/FJSVxtclanga/tcsds-latest/lib64:/opt/FJSVxtclanga/tcsds-latest/lib/:/vol0004/ra000007/data/jruiz/comp_libs_fujitsu/netcdf/lib:$LD_LIBRARY_PATH
  export OMP_NUM_THREADS=${OMP_NUM_THREADS:-48}
  export PARALLEL=$OMP_NUM_THREADS
  export OMP_STACKSIZE=${OMP_STACKSIZE:-2G}
fi

fail () { echo "ERROR: $*"; exit 1; }
[ -e "$OUT" ] && fail "$OUT exists"
[ -x "$LETKF_EXE" ] || fail "no $LETKF_EXE (./make.sh FUJITSU_FUGAKU)"
[ -x "$SYNTH_EXE" ] || fail "no $SYNTH_EXE (./make.sh FUJITSU_FUGAKU)"
[ $TRUTH -ge 1 -a $TRUTH -le $NMEM ] || fail "TRUTH must be in 1..$NMEM"
for m in $(seq 1 $NMEM); do
  [ -r $GUES/gues$(printf %05d $m) ] || fail "missing $GUES/gues$(printf %05d $m)"
done

#Ensemble members (all but the truth), numbered 1..NBV for the LETKF
ENS=$(seq 1 $NMEM | grep -vx $TRUTH)
echo "== OSSE analysis: truth member $TRUTH, $NBV ensemble members, ntemp=$NTEMP alpha_s=$ALPHA_S relax_mode=$RELAX"
echo "   guesses: $GUES"
echo "   letkf:   $LETKF_EXE ($(cd $ROOT && git describe --always --dirty 2>/dev/null))"

mkdir -p $OUT/work
command -v lfs > /dev/null && lfs setstripe -c 8 $OUT/work
cd $OUT/work
#Copies, not links: the LETKF writes the analyses over its input files.
i=0
for m in $ENS; do
  i=$((i+1))
  cp $GUES/gues$(printf %05d $m) gs01$(printf %05d $i) &
done
wait
ln -s $GUES/gues$(printf %05d $TRUTH) truth
cp gs01$(printf %05d $NBV) guesemean
cp gs01$(printf %05d $NBV) analemean

#Synthetic observations (src/osse/synth_obs.f90)
cat > synth_obs.namelist <<EOF
&synth_obs
  truth_file = 'truth'
  ens_prefix = 'gs01'
  nens = $NBV
  seed = $SEED
  obs_file = 'rad0101.dat'
  truth_obs_file = 'truth_rad0101.dat'
/
EOF
$MPIEXE $SYNTH_EXE > synth_obs.log 2>&1 || { tail -20 synth_obs.log; fail "synth_obs failed"; }
grep -E "Candidate points|Written|obshdxf" SYNO-00000 | sed 's/^/   /'

#LETKF
sed -e "s/__NBV__/$NBV/" -e "s/__NRADAR__/1/" -e "s/__NTEMP__/$NTEMP/" \
    -e "s/__ALPHA_S__/$ALPHA_S/" -e "s/__TEMP_RELAX_MODE__/$RELAX/" $NAMELIST > letkf.namelist
grep -q "__" letkf.namelist && fail "placeholders left in letkf.namelist"
$MPIEXE $LETKF_EXE > letkf.log 2>&1 || { tail -20 NOUT-00000; fail "letkf failed"; }
grep -q "LETKF ends" NOUT-00000 || { tail -20 NOUT-00000; fail "letkf did not reach the end"; }
grep -E "TEMPERING STEP|TIMER\(DAS_LETKF" NOUT-00000 | sed 's/^/   /'

#Results under the original member numbers
i=0
for m in $ENS; do
  i=$((i+1))
  mv gs01$(printf %05d $i) $OUT/anal$(printf %05d $m)
done
cp $GUES/gues$(printf %05d $TRUTH) $OUT/anal$(printf %05d $TRUTH)
mv guesemean analemean NOUT-00000 letkf.namelist rad0101.dat truth_rad0101.dat SYNO-00000 $OUT/
[ -e obs.dat ] && mv obs.dat $OUT/
cd $OUT
[ "${KEEP_WORK:-0}" == "1" ] || rm -rf $OUT/work
echo "== Done: $OUT"
ls $OUT | head -8
