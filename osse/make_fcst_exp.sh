#!/bin/bash
#=======================================================================
# Build a wrfutil DAFCST experiment that runs ensemble forecasts from a
# given set of analyses, reusing the scripts, configuration, namelists,
# WRF executables and boundary data (met_em) of an existing experiment.
# Nothing in the source experiment is modified: files are copied or linked.
#
# Usage:
#   osse/make_fcst_exp.sh [--dry-run] SRC_EXP NEW_EXP "YYYY-MM-DD hh:mm:ss" \
#                         LEAD_SECONDS MEM_INI MEM_END ANAL_DIR [INODE]
#
#   SRC_EXP   experiment to copy (e.g. the director's 2 km 5-min experiment)
#   NEW_EXP   new experiment directory (must not exist), in the group data
#             area: /data/ra000007/$USER/...
#   date      forecast initialization (UTC)
#   LEAD      forecast length [s]; outputs every 300 s
#   MEM_INI, MEM_END  members to run (3 digits, e.g. 001 060)
#   ANAL_DIR  directory with anal<NNNNN> files valid at the initialization;
#             linked as NEW_EXP/HIST/ANAL/<YYYYmmddHHMMSS>
#   INODE     nodes for the job (default 240: 30 members of 8 nodes at a
#             time, which keeps the job in the "small" resource group)
#
# OSSE jobs (thesis Chapter 4 extension):
#   job 1: init 2024-03-19 17:55:00, LEAD 3900, members 001-060, analyses of
#          the source experiment at 17:55 -> prior at 18:00, truth (060) and
#          forecast without assimilation up to 19:00.
#   job 3: init 2024-03-19 18:00:00, LEAD 3600, members 001-059, the LETKF
#          analyses at 18:00.
#
# Boundary data: EXTWPSPATH points to SRC_EXP/HIST/WPS/met_em, so main_FCST.sh
# links it at step 0 and nothing is copied. A lead time that is not a
# multiple of BDY_FREQ (1800 s) makes run_ensfcst.sh interpolate met_em in
# time (as in the 5-min DA cycle): LEAD 3900 does, LEAD 3600 does not.
#
# --dry-run checks every input path and prints the configuration changes
# without creating anything.
#
# After building: cd NEW_EXP/bin && ./queue_sub.sh main_FCST.sh
# Forecasts land in NEW_EXP/HIST/DAFCST/<YYYYmmddHHMMSS>/<MEM>/.
# Afterwards remove the work directories NEW_EXP/WRF/<MEM>/ (links,
# interpolated met_em, wrfbdy): the group is close to its inode limit.
#=======================================================================
set -e
DRY=0
if [ "$1" == "--dry-run" ]; then DRY=1; shift; fi
if [ $# -lt 7 ]; then sed -n '3,/^#====/p' "$0"; exit 1; fi

SRC=$(readlink -f "$1")
NEW=$2
INIT=$3
LEAD=$4
MEM_INI=$5
MEM_END=$6
ANAL_DIR=$(readlink -f "$7")
INODE=${8:-240}
INIT_TAG=$(date -u -d "$INIT UTC" +%Y%m%d%H%M%S)
MET_EM=$SRC/HIST/WPS/met_em
nfail=0

check () {   # $1 = description, $2... = test
  local d=$1; shift
  if "$@" ; then echo "  ok    $d"; else echo "  FAIL  $d"; nfail=$((nfail+1)); fi
}

echo "== Checks"
for d in bin lib conf NAMELIST WRF/code; do
  check "readable $SRC/$d" test -r $SRC/$d -a -x $SRC/$d
done
check "real.exe, wrf.exe, da_update_bc.exe in WRF/code" \
      test -e $SRC/WRF/code/real.exe -a -e $SRC/WRF/code/wrf.exe -a -e $SRC/WRF/code/da_update_bc.exe
check "main_FCST.sh handles EXTWPSPATH" grep -q EXTWPSPATH $SRC/bin/main_FCST.sh
check "queue_sub.sh present" test -x $SRC/bin/queue_sub.sh
check "readable met_em ($MET_EM)" test -r $MET_EM -a -x $MET_EM
nmiss=0
for m in $(seq $((10#$MEM_INI)) $((10#$MEM_END))); do
  f=$ANAL_DIR/anal$(printf %05d $m)
  [ -r $f ] || { echo "        missing $f"; nmiss=$((nmiss+1)); }
done
check "analyses anal$(printf %05d $((10#$MEM_INI)))..anal$(printf %05d $((10#$MEM_END))) in $ANAL_DIR" test $nmiss -eq 0
check "new experiment does not exist yet" test ! -e $NEW
case $(readlink -m $NEW) in
  /vol0004/*|/data/*) check "new experiment in the data area" true ;;
  *) [ -n "$ALLOW_ANY_DIR" ] && check "new experiment outside the data area (ALLOW_ANY_DIR, tests only)" true && true ;;&
  *) [ -z "$ALLOW_ANY_DIR" ] && check "new experiment in the data area (got $NEW)" false ;;
esac

#Boundary dates the forecast will need (same rules as lib/timeutils.sh)
BDY_FREQ=$(grep -E '^ *(export +)?BDY_FREQ=' $SRC/conf/exp.conf | head -1 | sed 's/.*BDY_FREQ=\([0-9]*\).*/\1/')
BDY_INI_FREQ=$(grep -E '^ *(export +)?BDY_INI_FREQ=' $SRC/conf/exp.conf | head -1 | sed 's/.*BDY_INI_FREQ=\([0-9]*\).*/\1/')
t0=$(date -u -d "$INIT UTC" +%s)
INI_BDY_TAG=$(date -u -d @$(( t0 / BDY_INI_FREQ * BDY_INI_FREQ )) +%Y%m%d%H%M%S)
echo "  info  BDY_FREQ=$BDY_FREQ s, BDY_INI_FREQ=$BDY_INI_FREQ s -> met_em/$INI_BDY_TAG/<MEM>/"
if (( LEAD % BDY_FREQ != 0 )); then
  echo "  info  LEAD=$LEAD is not a multiple of BDY_FREQ: met_em will be interpolated in time"
fi
if [ -r $MET_EM -a -x $MET_EM ]; then
  check "met_em/$INI_BDY_TAG/$MEM_INI exists" test -d $MET_EM/$INI_BDY_TAG/$MEM_INI
fi

echo "== Configuration changes"
cat <<EOF
  conf/config.env : BASEDIR=$(readlink -m $NEW)/
  conf/exp.conf   : EXPTYPE=DAFCST  RUN_WPS=0  RUN_BDY_PERT=0  RUN_LETKF=0  RUN_FCST=1
                    EXTWPSPATH=$MET_EM
                    BDY_MEM_INI=$MEM_INI  BDY_MEM_END=$MEM_END
                    FCST_INI_DATE=FCST_END_DATE='$INIT'  FCST_INI_FREQ=300
                    FCST_LEAD_TIME=$LEAD  FCST_OFREQ=300
  conf/machine.conf: INODE=$INODE  TOTAL_TIME_LIMIT=02:00:00
  HIST/ANAL/$INIT_TAG -> $ANAL_DIR
  WRF/code/*  -> links to $SRC/WRF/code/*
EOF

if [ $DRY -eq 1 ]; then
  echo "== Dry run: $nfail failed checks, nothing created"
  exit $nfail
fi
[ $nfail -eq 0 ] || { echo "== $nfail failed checks, nothing created"; exit 1; }

echo "== Building $NEW"
mkdir -p $NEW
NEW=$(readlink -f $NEW)
#Scripts and configuration, without job logs and caches
for d in bin lib conf NAMELIST; do
  cp -r $SRC/$d $NEW/
done
find $NEW/bin -maxdepth 1 \( -name '*.out' -o -name '*.err' -o -name '*.stats' -o -name '.stdin.*' \) -delete
find $NEW/bin -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
find $NEW/bin -name 'nohup.out' -delete
rm -f $NEW/conf/step.conf
mkdir -p $NEW/WRF/code $NEW/HIST/ANAL $NEW/HIST/DAFCST $NEW/LOGS $NEW/PROCS
ln -s $SRC/WRF/code/* $NEW/WRF/code/
ln -s $ANAL_DIR $NEW/HIST/ANAL/$INIT_TAG
if command -v lfs > /dev/null ; then lfs setstripe -c 8 $NEW/HIST/DAFCST ; fi

setvar () {  # $1 = file, $2 = variable, $3 = value (written verbatim)
  local f=$1 v=$2 val=$3
  if grep -qE "^ *(export +)?$v=" $f ; then
    sed -i -E "s|^ *(export +)?$v=.*|export $v=$val   # set by WRF_TEnKF osse/make_fcst_exp.sh|" $f
  else
    echo "export $v=$val   # set by WRF_TEnKF osse/make_fcst_exp.sh" >> $f
  fi
}
setvar $NEW/conf/config.env  BASEDIR         "$NEW/"
setvar $NEW/conf/exp.conf    EXPTYPE         '"DAFCST"'
setvar $NEW/conf/exp.conf    RUN_WPS         0
setvar $NEW/conf/exp.conf    RUN_BDY_PERT    0
setvar $NEW/conf/exp.conf    RUN_LETKF       0
setvar $NEW/conf/exp.conf    RUN_FCST        1
setvar $NEW/conf/exp.conf    EXTWPSPATH      "\"$MET_EM\""
setvar $NEW/conf/exp.conf    BDY_MEM_INI     $MEM_INI
setvar $NEW/conf/exp.conf    BDY_MEM_END     $MEM_END
setvar $NEW/conf/exp.conf    FCST_INI_DATE   "'$INIT'"
setvar $NEW/conf/exp.conf    FCST_END_DATE   "'$INIT'"
setvar $NEW/conf/exp.conf    FCST_INI_FREQ   300
setvar $NEW/conf/exp.conf    FCST_LEAD_TIME  $LEAD
setvar $NEW/conf/exp.conf    FCST_OFREQ      300
setvar $NEW/conf/machine.conf INODE          $INODE
setvar $NEW/conf/machine.conf TOTAL_TIME_LIMIT '"02:00:00"'

echo "== Changed lines"
grep -h "make_fcst_exp.sh" $NEW/conf/*.conf $NEW/conf/config.env | sed 's/   # set by.*//; s/^/  /'
echo "== Done. Submit with: cd $NEW/bin && ./queue_sub.sh main_FCST.sh"
