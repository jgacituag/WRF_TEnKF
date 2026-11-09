#!/bin/bash
#=======================================================================
# Build a wrfutil DA-cycle experiment (main_ASIM.sh) that repeats an
# existing experiment over a new period: same scripts, configuration,
# namelists, executables, boundary data (GEFS) and observations. Only the
# dates, BASEDIR and the experiment type change. Nothing in the source
# experiment is modified: files are copied or linked.
#
# Usage:
#   osse/make_da_exp.sh [--dry-run] [--smoke] [--time HH:MM:SS] SRC_EXP NEW_EXP \
#                       "DA_INI_DATE" "DA_END_DATE" [SPIN_UP_LENGTH]
#
#   SRC_EXP      experiment to copy (the director's 4 km 5-min experiment)
#   NEW_EXP      new experiment directory (must not exist), in the group
#                data area: /data/ra000007/$USER/...
#   DA_INI_DATE  start of the spin-up forecast (UTC, "YYYY-MM-DD hh:mm:ss")
#   DA_END_DATE  time of the last analysis (see "Steps" below)
#   SPIN_UP_LENGTH  [s], default: the value in SRC_EXP/conf/exp.conf
#
#   --smoke      4 members (BDY 001-004 -> target 001-004), INODE=8,
#                WRFNODE=2, NITER=1, TOTAL_TIME_LIMIT=00:30:00, interactive job:
#                a quick check of the whole chain in the "int" group.
#   --time       job time limit (default 24:00:00; 00:30:00 with --smoke)
#
# Steps (bin/main_ASIM.sh, lib/timeutils.sh): step 0 is the spin-up
# forecast DA_INI_DATE -> DA_INI_DATE + SPIN_UP_LENGTH, whose output is the
# first "analysis". Step s >= 1 is a guess from t_{s-1} to t_s = t_0 + s*300
# followed by the LETKF at t_s. The number of steps is
#   (DA_END_DATE - DA_INI_DATE - SPIN_UP_LENGTH)/ANALYSIS_FREQ + 1
# so the last analysis is at DA_END_DATE.
#
# Example (thesis OSSE, option B: 6 h spin-up from GEFS 06z, DA 12->18 UTC):
#   osse/make_da_exp.sh $SRC /data/ra000007/$USER/osse/da4km \
#        "2024-03-19 06:00:00" "2024-03-19 18:00:00"
#
# --dry-run checks every input (executables, GEFS files for every step,
# radar files) and prints the configuration changes without creating
# anything.
# After building: cd NEW_EXP/bin && ./queue_sub.sh main_ASIM.sh
# (or the job script of your choice that calls main_ASIM.sh).
# To restart after a time limit, set INI_STEP in NEW_EXP/conf/exp.conf to
# the first step that did not finish and submit again.
#=======================================================================
set -e
DRY=0 ; SMOKE=0 ; TLIM=""
while [[ "$1" == --* ]]; do
  case $1 in
    --dry-run) DRY=1 ; shift ;;
    --smoke)   SMOKE=1 ; shift ;;
    --time)    TLIM=$2 ; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done
if [ $# -lt 4 ]; then sed -n '3,/^#====/p' "$0"; exit 1; fi

SRC=$(readlink -f "$1")
NEW=$2
DA_INI=$3
DA_END=$4
NEWABS=$(readlink -m $NEW)
nfail=0

check () {   # $1 = description, $2... = test
  local d=$1; shift
  if "$@" ; then echo "  ok    $d"; else echo "  FAIL  $d"; nfail=$((nfail+1)); fi
}
confval () { # $1 = file, $2 = variable: first value assigned (quotes removed)
  grep -E "^ *(export +)?$2=" $1 | head -1 | sed -E "s/^ *(export +)?$2=//; s/[[:space:]]+#.*//; s/[[:space:]]+$//; s/^['\"]//; s/['\"]$//"
}

SPIN=${5:-$(confval $SRC/conf/exp.conf SPIN_UP_LENGTH)}
AFREQ=$(confval $SRC/conf/exp.conf ANALYSIS_FREQ)
AWIN=$(confval $SRC/conf/exp.conf ANALYSIS_WIN_END)
BDY_FREQ=$(confval $SRC/conf/exp.conf BDY_FREQ)
BDY_INI_FREQ=$(confval $SRC/conf/exp.conf BDY_INI_FREQ)
WPS_CYCLE=$(confval $SRC/conf/exp.conf WPS_CYCLE)
BDYPATH=$(confval $SRC/conf/exp.conf BDYPATH)
BDYPREFIX=$(confval $SRC/conf/exp.conf BDYPREFIX)
BDYVTABLE=$(confval $SRC/conf/exp.conf BDYVTABLE)
SRC_SOURCE=$(confval $SRC/conf/exp.conf WPS_DATA_SOURCE)
if [ $SMOKE -eq 1 ]; then
  BMI=001 ; BME=004 ; TMI=001 ; TME=004 ; INODE=8 ; WRFNODE=2 ; TLIM=${TLIM:-00:30:00} ; INTER=1
else
  BMI=$(confval $SRC/conf/exp.conf BDY_MEM_INI) ; BME=$(confval $SRC/conf/exp.conf BDY_MEM_END)
  TMI=$(confval $SRC/conf/exp.conf MEM_INI_TAR) ; TME=$(confval $SRC/conf/exp.conf MEM_END_TAR)
  INODE=$(confval $SRC/conf/machine.conf INODE) ; WRFNODE=$(confval $SRC/conf/machine.conf WRFNODE)
  TLIM=${TLIM:-24:00:00} ; INTER=0
fi

T_INI=$(date -u -d "$DA_INI UTC" +%s)
T_END=$(date -u -d "$DA_END UTC" +%s)
NSTEPS=$(( (T_END - T_INI - SPIN) / AFREQ + 1 ))
tag () { date -u -d @$1 +%Y%m%d%H%M%S; }

echo "== Period"
echo "  info  spin-up $(tag $T_INI) -> $(tag $((T_INI+SPIN))) ($SPIN s), $NSTEPS steps (0..$((NSTEPS-1))),"
echo "        last analysis $(tag $((T_INI+SPIN+(NSTEPS-1)*AFREQ)))"
echo "  info  members: boundary $BMI-$BME, target $TMI-$TME; INODE=$INODE WRFNODE=$WRFNODE limit $TLIM"

echo "== Checks"
check "SPIN_UP_LENGTH is a multiple of ANALYSIS_FREQ ($SPIN, $AFREQ)" test $((SPIN % AFREQ)) -eq 0
check "at least one DA step after the spin-up" test $NSTEPS -ge 2
check "source experiment takes boundaries from GFS/GEFS (WPS_DATA_SOURCE=$SRC_SOURCE)" test "$SRC_SOURCE" == "GFS"
check "WPS_CYCLE=1 in the source experiment" test "$WPS_CYCLE" == "1"
for d in bin lib conf NAMELIST WRF/code WPS/code PERT/code LETKF/code; do
  check "readable $SRC/$d" test -r $SRC/$d -a -x $SRC/$d
done
for f in WRF/code/real.exe WRF/code/wrf.exe WRF/code/da_update_bc.exe WPS/code/ungrib.exe \
         WPS/code/metgrid.exe PERT/code/pert_met_em.exe LETKF/code/letkf.exe \
         LETKF/code/update_wrf_time.exe WPS/geogrid/geo_em.d01.nc WPS/$BDYVTABLE; do
  check "present $f" test -e $SRC/$f
done
for f in bin/main_ASIM.sh bin/run_WPS.sh bin/run_Pert.sh bin/run_ensfcst.sh bin/run_LETKF.sh \
         bin/queue_sub.sh lib/timeutils.sh lib/encolarPJM.sh conf/obs.conf; do
  check "readable $f" test -r $SRC/$f
done
#Job logs, caches and hidden scratch files (some are unreadable by the group)
TAR_EXCLUDE="--exclude=*.out --exclude=*.err --exclude=*.stats --exclude=.* --exclude=__pycache__ --exclude=nohup.out"
UNREADABLE=$(cd $SRC && find bin lib conf NAMELIST \( -name '*.out' -o -name '*.err' -o -name '*.stats' \
              -o -name '.*' -o -name '__pycache__' -o -name 'nohup.out' \) -prune -o ! -readable -print 2>/dev/null)
[ -n "$UNREADABLE" ] && echo "$UNREADABLE" | sed 's/^/        unreadable /'
check "every file to copy from bin, lib, conf, NAMELIST is readable" test -z "$UNREADABLE"
check "new experiment does not exist yet" test ! -e $NEW
case $NEWABS in
  /vol0004/*|/data/*) check "new experiment in the data area" true ;;
  *) if [ -n "$ALLOW_ANY_DIR" ]; then check "new experiment outside the data area (ALLOW_ANY_DIR, tests only)" true
     else check "new experiment in the data area (got $NEW)" false ; fi ;;
esac

#GEFS files needed by every step (lib/timeutils.sh with WPS_CYCLE=1)
declare -A NEED
for (( s=0; s<NSTEPS; s++ )); do
  if [ $s -eq 0 ]; then t0=$T_INI ; t1=$((T_INI+SPIN))
  else t0=$((T_INI+SPIN+(s-1)*AFREQ)) ; t1=$((t0+AWIN)) ; fi
  ib=$(( t0 / BDY_INI_FREQ * BDY_INI_FREQ ))
  b0=$(( t0 / BDY_FREQ * BDY_FREQ )) ; b1=$(( (t1 + BDY_FREQ - 1) / BDY_FREQ * BDY_FREQ ))
  for (( b=b0; b<=b1; b+=BDY_FREQ )); do NEED["$ib $(( (b-ib)/3600 ))"]=1 ; done
done
nmiss=0 ; nfiles=0
for key in "${!NEED[@]}"; do
  read ib fh <<< "$key"
  dir=$BDYPATH/gefs.$(date -u -d @$ib +%Y%m%d)/$(date -u -d @$ib +%H)/$BDYPREFIX
  for m in $(seq -f %03g $((10#$BMI)) $((10#$BME))); do
    f=$(ls $dir/$m/g*pgrb2*.f$(printf %03d $fh) 2>/dev/null | head -1)
    if [ -n "$f" ] && [ -r "$f" ]; then nfiles=$((nfiles+1))
    else echo "        missing $dir/$m/g*pgrb2*.f$(printf %03d $fh)"; nmiss=$((nmiss+1)); fi
  done
done
echo "  info  GEFS (init, lead h): $(for k in "${!NEED[@]}"; do read ib fh <<< "$k"; echo "$(date -u -d @$ib +%d%H)z+$fh"; done | sort | tr '\n' ' ')"
check "GEFS files for every step ($nfiles found)" test $nmiss -eq 0

#Radar observations (informative: the LETKF runs without them)
OBSPATH=$(confval $SRC/conf/obs.conf OBSPATH)
RADARS=$(bash -c "source $SRC/conf/obs.conf > /dev/null 2>&1; echo \${RADARC_LIST[@]}")
nwith=0
for (( s=1; s<NSTEPS; s++ )); do
  ct=$(tag $((T_INI+SPIN+s*AFREQ)))
  ls $OBSPATH/RADARC/WRF_RADARC_*_$ct.dat > /dev/null 2>&1 && nwith=$((nwith+1))
done
echo "  info  radar files (OBSPATH=$OBSPATH, radars: $RADARS): $nwith of $((NSTEPS-1)) analysis times have at least one"

echo "== Configuration changes"
cat <<EOF
  conf/config.env  : BASEDIR=$NEWABS/
  conf/exp.conf    : EXPTYPE="DACYCLE"  INI_STEP=0  RUN_WPS=1  EXTWPSPATH=""  SAVEGUESS=1
                     RUN_FCST=1  RUN_LETKF=1  DA_INI_DATE='$DA_INI'  DA_END_DATE='$DA_END'
                     SPIN_UP_LENGTH=$SPIN  BDY_MEM_INI=$BMI  BDY_MEM_END=$BME
                     MEM_INI_TAR=$TMI  MEM_END_TAR=$TME
  conf/machine.conf: INODE=$INODE  WRFNODE=$WRFNODE  TOTAL_TIME_LIMIT="$TLIM"  INTERACTIVE_JOB=$INTER
  WRF/code, WPS/code, PERT/code, LETKF/code: links to the files of $SRC (namelists excluded)
  WPS/geogrid/geo_em.d01.nc -> $SRC/WPS/geogrid/
  WPS/Vtable* -> $SRC/WPS/ (ungrib reads \$WPSDIR/\$BDYVTABLE = $BDYVTABLE)
EOF

if [ $DRY -eq 1 ]; then
  echo "== Dry run: $nfail failed checks, nothing created"
  exit $nfail
fi
[ $nfail -eq 0 ] || { echo "== $nfail failed checks, nothing created"; exit 1; }

echo "== Building $NEWABS"
mkdir -p $NEWABS
NEW=$NEWABS
( cd $SRC && tar -cf - $TAR_EXCLUDE bin lib conf NAMELIST ) | ( cd $NEW && tar -xf - )
rm -f $NEW/conf/step.conf
mkdir -p $NEW/HIST/ANAL $NEW/HIST/GUES $NEW/HIST/WPS $NEW/LOGS $NEW/PROCS $NEW/WPS/geogrid
#Executables: link file by file, but never the namelists, which the run
#scripts overwrite with cp (cp would write through the link).
for d in WRF/code WPS/code PERT/code LETKF/code; do
  mkdir -p $NEW/$d
  for f in $SRC/$d/*; do
    case $(basename $f) in *namelist*) continue ;; esac
    ln -s $f $NEW/$d/
  done
done
ln -s $SRC/WPS/geogrid/geo_em.d01.nc $NEW/WPS/geogrid/
ln -s $SRC/WPS/Vtable* $NEW/WPS/          #run_WPS.sh links \$WPSDIR/\$BDYVTABLE as the ungrib Vtable
stripe () { if command -v lfs > /dev/null ; then lfs setstripe -c ${2:-8} $1 ; fi ; }
stripe $NEW/HIST/ANAL ; stripe $NEW/HIST/GUES ; stripe $NEW/HIST/WPS 4 ; mkdir -p $NEW/WRF && stripe $NEW/WRF 4

setvar () {  # $1 = file, $2 = variable, $3 = value (written verbatim)
  local f=$1 v=$2 val=$3
  if grep -qE "^ *(export +)?$v=" $f ; then
    sed -i -E "s|^ *(export +)?$v=.*|export $v=$val   # set by WRF_TEnKF osse/make_da_exp.sh|" $f
  else
    echo "export $v=$val   # set by WRF_TEnKF osse/make_da_exp.sh" >> $f
  fi
}
setvar $NEW/conf/config.env   BASEDIR         "$NEW/"
setvar $NEW/conf/exp.conf     EXPTYPE         '"DACYCLE"'
setvar $NEW/conf/exp.conf     INI_STEP        0
setvar $NEW/conf/exp.conf     RUN_WPS         1
setvar $NEW/conf/exp.conf     EXTWPSPATH      '""'
setvar $NEW/conf/exp.conf     SAVEGUESS       1
setvar $NEW/conf/exp.conf     RUN_FCST        1
setvar $NEW/conf/exp.conf     RUN_LETKF       1
setvar $NEW/conf/exp.conf     DA_INI_DATE     "'$DA_INI'"
setvar $NEW/conf/exp.conf     DA_END_DATE     "'$DA_END'"
setvar $NEW/conf/exp.conf     SPIN_UP_LENGTH  $SPIN
setvar $NEW/conf/exp.conf     BDY_MEM_INI     $BMI
setvar $NEW/conf/exp.conf     BDY_MEM_END     $BME
setvar $NEW/conf/exp.conf     MEM_INI_TAR     $TMI
setvar $NEW/conf/exp.conf     MEM_END_TAR     $TME
setvar $NEW/conf/machine.conf INODE           $INODE
setvar $NEW/conf/machine.conf WRFNODE         $WRFNODE
setvar $NEW/conf/machine.conf TOTAL_TIME_LIMIT "\"$TLIM\""
setvar $NEW/conf/machine.conf INTERACTIVE_JOB $INTER
#pert_met_em computes the members in NITER batches of ceil(nbv_tar/NITER)
#and runs past nbv_tar when NITER > nbv_tar (10 in the source, 4 members
#in the smoke test): a single batch is the same computation.
[ $SMOKE -eq 1 ] && setvar $NEW/conf/exp.conf NITER 1

echo "== Changed lines"
grep -h "make_da_exp.sh" $NEW/conf/*.conf $NEW/conf/config.env | sed 's/   # set by.*//; s/^/  /'
echo "== Done. Submit with: cd $NEW/bin && ./queue_sub.sh main_ASIM.sh"
