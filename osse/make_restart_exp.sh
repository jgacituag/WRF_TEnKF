#!/bin/bash
#=======================================================================
# Build a wrfutil DA-cycle experiment that continues an existing cycle from
# one of its analyses with the tempering LETKF of this repository. The new
# experiment starts at the step after RESTART_DATE, from the source's
# analyses at that time, and runs to the source's DA_END_DATE with the same
# scripts, configuration, observations, WRF executables and perturbed
# boundaries (the source's HIST/WPS/met_em, so WPS and pert_met_em do not
# run). Only the LETKF changes: this repository's executables (a copy of
# the tar), bin/run_LETKF.sh and namelists/letkf.namelist, which add the
# &tempering block, with NTEMP, ALPHA_S and TEMP_RELAX_MODE in letkf.conf.
# The source's inflation, relaxation, localization and observations stay.
# With --ntemp 1 the LETKF reproduces the source's bit for bit (regression
# test), so that run repeats the source.
#
# Nothing in the source is modified: its analyses at RESTART_DATE and its
# met_em are linked file by file (cp -rs) into folders of the new
# experiment, so no step can write into the source.
#
# Usage:
#   osse/make_restart_exp.sh [--dry-run] [options] SRC_EXP NEW_EXP "RESTART_DATE"
#
#   SRC_EXP       a DA cycle with HIST/ANAL/<RESTART_DATE> complete and
#                 HIST/WPS/met_em (e.g. the 4 km, 5-min run A)
#   NEW_EXP       new experiment directory (must not exist)
#   RESTART_DATE  analysis time to start from (UTC, "YYYY-MM-DD hh:mm:ss");
#                 the first new step forecasts from it to the next analysis
#
# Options:
#   --ntemp N         tempering steps [2]
#   --alpha-s A       delay parameter of the weights [1.0d0]
#   --relax-mode M    TEMP_RELAX_MODE, 1 or 2 [1]
#   --letkf TAR       LETKF executables [dist/letkf_<BUILD>.tar of this repository]
#   --time HH:MM:SS   job time limit [the source's TOTAL_TIME_LIMIT]
#
# --dry-run checks every input (the analyses, the met_em of every member,
# observations at every step, the LETKF tar, that this repository's
# run_LETKF.sh and letkf.namelist differ from the source's only by the
# tempering lines) and prints the changes without creating anything.
#
# The new experiment keeps every guess and analysis from RESTART_DATE on
# (keep_times.txt, read by prune_hist.sh). Submit with the source's job
# script, e.g.  cd NEW_EXP/bin && ./queue_sub.sh main_ASIM_keep_retry.sh
# After the job: NEW_EXP/clean_work.sh removes the work folders.
#
# Example (experiment E, real-data cycles from 17:00 with Nt = 2):
#   osse/make_restart_exp.sh --dry-run --ntemp 2 ~/salidas/osse/da4km \
#        ~/salidas/osse/e17/e_nt2 "2024-03-19 17:00:00"
#=======================================================================
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
DRY=0 ; NTEMP=2 ; ALPHA_S=1.0d0 ; RELAX=1 ; TLIM=""
BUILD=${BUILD:-$( [ -d /opt/FJSVxtclanga ] && echo FUJITSU_FUGAKU || echo INTEL_HYDRA )}
TAR=$ROOT/dist/letkf_$BUILD.tar
while [[ "$1" == --* ]]; do
  case $1 in
    --dry-run)    DRY=1 ; shift ;;
    --ntemp)      NTEMP=$2 ; shift 2 ;;
    --alpha-s)    ALPHA_S=$2 ; shift 2 ;;
    --relax-mode) RELAX=$2 ; shift 2 ;;
    --letkf)      TAR=$(readlink -f "$2") ; shift 2 ;;
    --time)       TLIM=$2 ; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done
if [ $# -ne 3 ]; then sed -n '3,/^#====/p' "$0"; exit 1; fi

SRC=$(readlink -f "$1")
NEW=$(readlink -m "$2")
RESTART=$3
TAG=$(date -u -d "$RESTART UTC" +%Y%m%d%H%M%S)
nfail=0

check () {   # $1 = description, $2... = test
  local d=$1; shift
  if "$@" ; then echo "  ok    $d"; else echo "  FAIL  $d"; nfail=$((nfail+1)); fi
}
confstr () { # $1 = file, $2 = variable: first value assigned, quotes and comment removed
  grep -E "^ *(export +)?$2=" $1 | head -1 | sed -E "s/^ *(export +)?$2=//; s/[[:space:]]+#.*//; s/[[:space:]]+$//; s/^['\"]//; s/['\"]$//"
}
tempering_only () {  # $1 = this repository's file, $2 = the source's: same but for the tempering lines
  diff -q <(grep -v -e '__NTEMP__' -e '__ALPHA_S__' -e '__TEMP_RELAX_MODE__' -e '^&tempering' \
                   -e 'Tempering (WRF_TEnKF)' $1 | sed '/^\/$/{N;s/^\/\n\//\//}') $2 > /dev/null 2>&1
}

#Source settings
DA_INI=$(confstr $SRC/conf/exp.conf DA_INI_DATE)
DA_END=$(confstr $SRC/conf/exp.conf DA_END_DATE)
SPIN=$(confstr $SRC/conf/exp.conf SPIN_UP_LENGTH)
FREQ=$(confstr $SRC/conf/exp.conf ANALYSIS_FREQ)
BDYFREQ=$(confstr $SRC/conf/exp.conf BDY_FREQ)
BDYINI=$(confstr $SRC/conf/exp.conf BDY_INI_FREQ)
MEM_INI=$(confstr $SRC/conf/exp.conf MEM_INI_TAR)
MEM_END=$(confstr $SRC/conf/exp.conf MEM_END_TAR)
OBSPATH=$(confstr $SRC/conf/obs.conf OBSPATH)
TLIM=${TLIM:-$(confstr $SRC/conf/machine.conf TOTAL_TIME_LIMIT)}
NMEM=$(( 10#$MEM_END - 10#$MEM_INI + 1 ))
T0=$(date -u -d "$RESTART UTC" +%s)
T1=$(date -u -d "$DA_END UTC" +%s)
S0=$(( (T0 - $(date -u -d "$DA_INI UTC" +%s) - SPIN) / FREQ ))
INI_STEP=$(( S0 + 1 ))
NSTEP=$(( (T1 - T0) / FREQ ))
STEP_TAGS=$(for ((t = T0 + FREQ; t <= T1; t += FREQ)); do date -u -d @$t +%Y%m%d%H%M%S; done)
INIB=$(( T0 / BDYINI * BDYINI ))                          #boundary initialization of the window
INIB_TAG=$(date -u -d @$INIB +%Y%m%d%H%M%S)
BFIRST=$(date -u -d @$(( T0 / BDYFREQ * BDYFREQ )) +%Y-%m-%d_%H:%M:%S)
BLAST=$(date -u -d @$(( (T1 + BDYFREQ - 1) / BDYFREQ * BDYFREQ )) +%Y-%m-%d_%H:%M:%S)

echo "== Checks"
check "new experiment does not exist" test ! -e $NEW
for d in bin lib conf NAMELIST WRF/code; do
  check "readable $SRC/$d" test -r $SRC/$d -a -x $SRC/$d
done
for f in bin/main_ASIM.sh bin/run_ensfcst.sh bin/run_LETKF.sh bin/main_ASIM_keep_retry.sh \
         bin/prune_hist.sh bin/clean_wrf_members.sh lib/timeutils.sh NAMELIST/letkf.namelist; do
  check "readable $f" test -r $SRC/$f
done
check "RESTART_DATE on the analysis grid (step $S0 of the source)" \
      test $(( (T0 - $(date -u -d "$DA_INI UTC" +%s) - SPIN) % FREQ )) -eq 0 -a $S0 -ge 1
check "RESTART_DATE before DA_END_DATE ($NSTEP steps to run)" test $NSTEP -ge 1
NA=$(ls $SRC/HIST/ANAL/$TAG/anal[0-9]* 2>/dev/null | wc -l)
check "$NMEM analyses in HIST/ANAL/$TAG (found $NA)" test $NA -eq $NMEM
check "HIST/ANAL/$TAG/NOUT-00000 says LETKF ends" grep -qs "LETKF ends" $SRC/HIST/ANAL/$TAG/NOUT-00000
NM=0
for m in $(seq -f %03g $((10#$MEM_INI)) $((10#$MEM_END))); do
  D=$SRC/HIST/WPS/met_em/$INIB_TAG/$m
  [ -r $D/met_em.d01.$BFIRST.nc -a -r $D/met_em.d01.$BLAST.nc ] && NM=$((NM+1))
done
check "met_em $BFIRST and $BLAST in HIST/WPS/met_em/$INIB_TAG for $NM/$NMEM members" test $NM -eq $NMEM
NOBS=0
for t in $STEP_TAGS; do ls $OBSPATH/*/WRF_*_$t.dat > /dev/null 2>&1 && NOBS=$((NOBS+1)); done
check "observation files in $OBSPATH at $NOBS/$NSTEP analysis times" test $NOBS -eq $NSTEP
check "LETKF tar $TAR holds letkf.exe" bash -c "tar -tf $TAR 2>/dev/null | grep -qx letkf.exe"
check "bin/run_LETKF.sh = the source's + tempering lines" tempering_only $ROOT/bin/run_LETKF.sh $SRC/bin/run_LETKF.sh
check "namelists/letkf.namelist = the source's + &tempering block" tempering_only $ROOT/namelists/letkf.namelist $SRC/NAMELIST/letkf.namelist
check "NTEMP=$NTEMP between 1 and 20, TEMP_RELAX_MODE=$RELAX 1 or 2" \
      test "$NTEMP" -ge 1 -a "$NTEMP" -le 20 -a \( "$RELAX" == 1 -o "$RELAX" == 2 \)

echo "== Configuration changes"
cat <<EOF
  conf/config.env  : BASEDIR=$NEW/  LETKFPATH=$NEW/LETKF/letkf.tar
  conf/exp.conf    : INI_STEP=$INI_STEP (first analysis $(echo "$STEP_TAGS" | head -1), last $(echo "$STEP_TAGS" | tail -1))
                     RUN_WPS=0  RUN_BDY_PERT=0  EXTWPSPATH=""
  conf/letkf.conf  : NTEMP=$NTEMP  ALPHA_S=$ALPHA_S  TEMP_RELAX_MODE=$RELAX
  conf/machine.conf: TOTAL_TIME_LIMIT=$TLIM
  bin/run_LETKF.sh, NAMELIST/letkf.namelist <- $ROOT (copies)
  LETKF/letkf.tar  <- $TAR (copy)
  HIST/ANAL/$TAG/*  -> $SRC/HIST/ANAL/$TAG/* (file links)
  HIST/WPS/met_em/*/*/* -> $SRC/HIST/WPS/met_em (file links)
  WRF/code/*       -> $SRC/WRF/code/*
  keep_times.txt   : $TAG and the $NSTEP new analysis times
EOF

if [ $DRY -eq 1 ]; then
  echo "== Dry run: $nfail failed checks, nothing created"
  exit $nfail
fi
[ $nfail -eq 0 ] || { echo "== $nfail failed checks, nothing created"; exit 1; }

echo "== Building $NEW"
mkdir -p $NEW
TAR_EXCLUDE="--exclude=*.out --exclude=*.err --exclude=*.stats --exclude=.* --exclude=__pycache__
             --exclude=nohup.out --exclude=*.o[0-9]* --exclude=*.e[0-9]* --exclude=bin/old"
( cd $SRC && tar -cf - $TAR_EXCLUDE bin lib conf NAMELIST ) | ( cd $NEW && tar -xf - )
rm -f $NEW/conf/step.conf
cp $ROOT/bin/run_LETKF.sh $NEW/bin/run_LETKF.sh
cp $ROOT/namelists/letkf.namelist $NEW/NAMELIST/letkf.namelist
mkdir -p $NEW/LETKF $NEW/WRF/code $NEW/HIST/ANAL $NEW/HIST/WPS $NEW/LOGS $NEW/PROCS
cp $TAR $NEW/LETKF/letkf.tar              #a copy, so a later rebuild cannot change a running experiment
ln -s $SRC/WRF/code/* $NEW/WRF/code/
cp -rs $SRC/HIST/ANAL/$TAG $NEW/HIST/ANAL/$TAG
cp -rs $SRC/HIST/WPS/met_em $NEW/HIST/WPS/met_em
{ echo $TAG ; echo "$STEP_TAGS" ; } > $NEW/keep_times.txt

setvar () {  # $1 = file, $2 = variable, $3 = value (written verbatim)
  local f=$1 v=$2 val=$3
  if grep -qE "^ *(export +)?$v=" $f ; then
    sed -i -E "s|^ *(export +)?$v=.*|export $v=$val   # set by WRF_TEnKF osse/make_restart_exp.sh|" $f
  else
    echo "export $v=$val   # set by WRF_TEnKF osse/make_restart_exp.sh" >> $f
  fi
}
setvar $NEW/conf/config.env   BASEDIR          "$NEW/"
setvar $NEW/conf/config.env   LETKFPATH        "$NEW/LETKF/letkf.tar"
setvar $NEW/conf/exp.conf     INI_STEP         $INI_STEP
setvar $NEW/conf/exp.conf     RUN_WPS          0
setvar $NEW/conf/exp.conf     RUN_BDY_PERT     0
setvar $NEW/conf/exp.conf     EXTWPSPATH       '""'
setvar $NEW/conf/letkf.conf   NTEMP            $NTEMP
setvar $NEW/conf/letkf.conf   ALPHA_S          $ALPHA_S
setvar $NEW/conf/letkf.conf   TEMP_RELAX_MODE  $RELAX
setvar $NEW/conf/machine.conf TOTAL_TIME_LIMIT "\"$TLIM\""

#Cleanup helper for after the job: work directories only, never outputs.
cat > $NEW/clean_work.sh <<EOF
#!/bin/bash
# Remove the work directories of this experiment (made by make_restart_exp.sh).
# Keeps HIST/ (guesses, analyses) and LOGS/.
cd $NEW || exit 1
rm -rf WRF/[0-9][0-9][0-9] LETKF/00 PROCS/*
echo "Removed work directories under $NEW"
EOF
chmod +x $NEW/clean_work.sh

echo "== Changed lines"
grep -h "make_restart_exp.sh" $NEW/conf/*.conf $NEW/conf/config.env | sed 's/   # set by.*//; s/^/  /'
echo "== Done. Submit with: cd $NEW/bin && ./queue_sub.sh main_ASIM_keep_retry.sh"
echo "   After the job: $NEW/clean_work.sh"
