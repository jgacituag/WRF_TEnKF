#!/bin/bash
#=======================================================================
# Build a wrfutil DAFCST experiment that runs ensemble forecasts from a
# given set of analyses, reusing the scripts, configuration, namelists and
# WRF/WPS executables of an existing experiment. Nothing in the source
# experiment is modified: files are copied or linked.
#
# Usage:
#   osse/make_fcst_exp.sh [--dry-run] [BOUNDARY OPTION] SRC_EXP NEW_EXP \
#                         "YYYY-MM-DD hh:mm:ss" LEAD MEM_INI MEM_END ANAL_DIR [INODE]
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
#   --time HH:MM:SS  job time limit (default 02:00:00)
#
# Boundary options (lateral boundaries come from met_em files):
#   (none)                    use SRC_EXP/HIST/WPS/met_em through EXTWPSPATH.
#   --gefs                    as the source DA experiment did: WPS from GEFS
#                             and pert_met_em inside the job (RUN_WPS=1,
#                             BDY_PERT and the boundary/target members of
#                             SRC_EXP, which must be a GFS experiment such as
#                             one built by make_da_exp.sh). pert_met_em needs
#                             the whole target ensemble, so MEM_INI..MEM_END
#                             must be the target members of SRC_EXP (for the
#                             OSSE, member 060 is the truth: its "analysis" is
#                             its guess). Same random seed and inputs as the DA
#                             run, so the same perturbed boundaries.
#   --bdy-from-anal ANAL_ROOT run WPS inside the job on wrfouts linked from
#                             ANAL_ROOT/<YYYYmmddHHMMSS>/anal<MEM> (each member
#                             takes its boundaries from its own analyses):
#                             WPS_DATA_SOURCE='WRF', wrf_to_wps + metgrid, with
#                             SRC_EXP/WPS/code and its geo_em. The met_em go to
#                             NEW_EXP/HIST/WPS/met_em_ori/. wrf_to_wps.exe is
#                             the WRF_TEnKF build, which adds the model-top
#                             level that real.exe needs (src/wps/wrf_to_wps.f90);
#                             METLEV grows by one accordingly. Override the
#                             executable with --wrf-to-wps EXE.
#   --met-em-links MET_ROOT   link existing met_em.d01.<date>.nc files found
#                             under MET_ROOT/*/<MEM>/ into the folder this
#                             forecast expects (e.g. job 3 reusing the met_em
#                             that job 1 produced).
#
# OSSE jobs (thesis Chapter 4 extension):
#   job 1: --bdy-from-anal <src>/HIST/ANAL, init 2024-03-19 17:55:00,
#          LEAD 3900, members 001-060, ANAL_DIR <src>/HIST/ANAL/20240319175500
#          -> prior at 18:00, truth (060) and forecast without DA up to 19:00.
#   job 3: --met-em-links <job1>/HIST/WPS/met_em_ori, init 2024-03-19 18:00:00,
#          LEAD 3600, members 001-059, ANAL_DIR = the LETKF analyses.
#
# Boundary times follow lib/timeutils.sh: from floor(init, BDY_FREQ) to
# ceil(init + LEAD, BDY_FREQ) every BDY_FREQ (1800 s), in the folder
# floor(init, BDY_INI_FREQ). A lead time or initialization off the BDY_FREQ
# grid makes run_ensfcst.sh interpolate met_em in time, as in the 5-min DA
# cycle.
#
# --dry-run checks every input and prints the configuration changes
# without creating anything.
#
# After building: cd NEW_EXP/bin && ./queue_sub.sh main_FCST.sh
# Forecasts land in NEW_EXP/HIST/DAFCST/<YYYYmmddHHMMSS>/<MEM>/.
# Afterwards run NEW_EXP/clean_work.sh to remove the work directories
# (links, intermediate files, interpolated met_em, wrfinput, wrfbdy): the
# group is close to its inode limit.
#=======================================================================
set -e
DRY=0 ; BDYMODE=src ; BDYARG="" ; TLIM=02:00:00
W2W=$(cd "$(dirname "$0")/.." && pwd)/build/FUJITSU_FUGAKU/wrf_to_wps.exe
while [[ "$1" == --* ]]; do
  case $1 in
    --dry-run)       DRY=1 ; shift ;;
    --bdy-from-anal) BDYMODE=anal ; BDYARG=$(readlink -f "$2") ; shift 2 ;;
    --met-em-links)  BDYMODE=links ; BDYARG=$(readlink -f "$2") ; shift 2 ;;
    --gefs)          BDYMODE=gefs ; shift ;;
    --time)          TLIM=$2 ; shift 2 ;;
    --wrf-to-wps)    W2W=$(readlink -f "$2") ; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done
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
NEWABS=$(readlink -m $NEW)
MEMS=$(seq -f %03g $((10#$MEM_INI)) $((10#$MEM_END)))
nfail=0

check () {   # $1 = description, $2... = test
  local d=$1; shift
  if "$@" ; then echo "  ok    $d"; else echo "  FAIL  $d"; nfail=$((nfail+1)); fi
}
confval () { # $1 = file, $2 = variable: first numeric value assigned
  grep -E "^ *(export +)?$2=" $1 | head -1 | sed -E "s/.*$2=['\"]?([0-9]*).*/\1/"
}
confstr () { # $1 = file, $2 = variable: first value assigned, quotes and comment removed
  grep -E "^ *(export +)?$2=" $1 | head -1 | sed -E "s/^ *(export +)?$2=//; s/[[:space:]]+#.*//; s/[[:space:]]+$//; s/^['\"]//; s/['\"]$//"
}

TAR_EXCLUDE="--exclude=*.out --exclude=*.err --exclude=*.stats --exclude=.* --exclude=__pycache__ --exclude=nohup.out"
echo "== Checks"
for d in bin lib conf NAMELIST WRF/code; do
  check "readable $SRC/$d" test -r $SRC/$d -a -x $SRC/$d
done
check "real.exe, wrf.exe, da_update_bc.exe in WRF/code" \
      test -e $SRC/WRF/code/real.exe -a -e $SRC/WRF/code/wrf.exe -a -e $SRC/WRF/code/da_update_bc.exe
check "main_FCST.sh handles EXTWPSPATH" grep -q EXTWPSPATH $SRC/bin/main_FCST.sh
UNREADABLE=$(cd $SRC && find bin lib conf NAMELIST \( -name '*.out' -o -name '*.err' -o -name '*.stats' \
              -o -name '.*' -o -name '__pycache__' -o -name 'nohup.out' \) -prune -o ! -readable -print 2>/dev/null)
[ -n "$UNREADABLE" ] && echo "$UNREADABLE" | sed 's/^/        unreadable /'
check "every file to copy from bin, lib, conf, NAMELIST is readable" test -z "$UNREADABLE"
for f in bin/main_FCST.sh bin/run_ensfcst.sh bin/run_WPS.sh lib/timeutils.sh lib/encolarPJM.sh; do
  check "readable $f" test -r $SRC/$f
done
check "queue_sub.sh present" test -x $SRC/bin/queue_sub.sh
nmiss=0
for m in $MEMS; do
  [ -r $ANAL_DIR/anal00$m ] || { echo "        missing $ANAL_DIR/anal00$m"; nmiss=$((nmiss+1)); }
done
check "initial analyses anal00$MEM_INI..anal00$MEM_END in $ANAL_DIR" test $nmiss -eq 0
check "new experiment does not exist yet" test ! -e $NEW
case $NEWABS in
  /vol0004/*|/data/*) check "new experiment in the data area" true ;;
  *) if [ -n "$ALLOW_ANY_DIR" ]; then check "new experiment outside the data area (ALLOW_ANY_DIR, tests only)" true
     else check "new experiment in the data area (got $NEW)" false ; fi ;;
esac

#Boundary times (lib/timeutils.sh)
BDY_FREQ=$(confval $SRC/conf/exp.conf BDY_FREQ)
BDY_INI_FREQ=$(confval $SRC/conf/exp.conf BDY_INI_FREQ)
WPS_CYCLE=$(confval $SRC/conf/exp.conf WPS_CYCLE)
t0=$(date -u -d "$INIT UTC" +%s) ; t1=$(( t0 + LEAD ))
BDY_T0=$(( t0 / BDY_FREQ * BDY_FREQ ))
BDY_T1=$(( (t1 + BDY_FREQ - 1) / BDY_FREQ * BDY_FREQ ))
if [ "$WPS_CYCLE" == "1" ]; then
  INI_BDY_TAG=$(date -u -d @$(( t0 / BDY_INI_FREQ * BDY_INI_FREQ )) +%Y%m%d%H%M%S)
else
  DA0=$(date -u -d "$(grep -E '^ *(export +)?DA_INI_DATE=' $SRC/conf/exp.conf | head -1 | cut -d"'" -f2) UTC" +%s)
  INI_BDY_TAG=$(date -u -d @$(( DA0 / BDY_INI_FREQ * BDY_INI_FREQ )) +%Y%m%d%H%M%S)
fi
BDY_DATES=$(for (( t=BDY_T0; t<=BDY_T1; t+=BDY_FREQ )); do date -u -d @$t +%Y%m%d%H%M%S; done)
echo "  info  BDY_FREQ=$BDY_FREQ s, folder met_em/$INI_BDY_TAG/<MEM>/, times: $(echo $BDY_DATES | sed 's/\([0-9]\{8\}\)\([0-9]\{4\}\)[0-9]\{2\}/\2/g')"
(( LEAD % BDY_FREQ != 0 || t0 % BDY_FREQ != 0 )) && \
  echo "  info  initialization or lead off the $BDY_FREQ s grid: met_em will be interpolated in time"
wpsname () { date -u -d "${1:0:8} ${1:8:2}:${1:10:2}:${1:12:2} UTC" +%Y-%m-%d_%H:%M:%S; }

case $BDYMODE in
src)
  MET_EM=$SRC/HIST/WPS/met_em
  check "readable met_em ($MET_EM)" test -r $MET_EM -a -x $MET_EM
  [ -r $MET_EM -a -x $MET_EM ] && check "met_em/$INI_BDY_TAG/$MEM_INI exists" test -d $MET_EM/$INI_BDY_TAG/$MEM_INI
  ;;
anal)
  check "readable WPS executables ($SRC/WPS/code)" test -r $SRC/WPS/code -a -x $SRC/WPS/code
  check "metgrid.exe in $SRC/WPS/code" test -e $SRC/WPS/code/metgrid.exe
  check "WRF_TEnKF wrf_to_wps.exe with model-top level ($W2W)" test -x $W2W
  METLEV_SRC=$(confval $SRC/conf/model.conf METLEV)
  echo "  info  METLEV $METLEV_SRC -> $((METLEV_SRC+1)) (extra model-top level)"
  check "METGRID.TBL in $SRC/WPS/code/metgrid" test -e $SRC/WPS/code/metgrid/METGRID.TBL
  check "readable geo_em ($SRC/WPS/geogrid/geo_em.d01.nc)" test -r $SRC/WPS/geogrid/geo_em.d01.nc
  check "main_FCST.sh links met_em_ori when BDY_PERT=0" grep -q met_em_ori $SRC/bin/main_FCST.sh
  nmiss=0
  for d in $BDY_DATES; do for m in $MEMS; do
    [ -r $BDYARG/$d/anal00$m ] || { echo "        missing $BDYARG/$d/anal00$m"; nmiss=$((nmiss+1)); }
  done; done
  check "boundary analyses ($(echo $BDY_DATES | wc -w) times x $(echo $MEMS | wc -w) members) in $BDYARG" test $nmiss -eq 0
  ;;
gefs)
  check "source experiment takes boundaries from GFS/GEFS" test "$(confstr $SRC/conf/exp.conf WPS_DATA_SOURCE)" == "GFS"
  TMI=$(confstr $SRC/conf/exp.conf MEM_INI_TAR) ; TME=$(confstr $SRC/conf/exp.conf MEM_END_TAR)
  BMI=$(confstr $SRC/conf/exp.conf BDY_MEM_INI) ; BME=$(confstr $SRC/conf/exp.conf BDY_MEM_END)
  check "forecast members $MEM_INI-$MEM_END are the target members of the source ($TMI-$TME)" \
        test "$MEM_INI" == "$TMI" -a "$MEM_END" == "$TME"
  for f in WPS/code/ungrib.exe WPS/code/metgrid.exe PERT/code/pert_met_em.exe WPS/geogrid/geo_em.d01.nc \
           WPS/$(confstr $SRC/conf/exp.conf BDYVTABLE); do
    check "present $f" test -e $SRC/$f
  done
  BDYPATH=$(confstr $SRC/conf/exp.conf BDYPATH) ; BDYPREFIX=$(confstr $SRC/conf/exp.conf BDYPREFIX)
  ib=$(date -u -d "${INI_BDY_TAG:0:8} ${INI_BDY_TAG:8:2}:00 UTC" +%s)
  nmiss=0 ; nfiles=0
  for d in $BDY_DATES; do
    td=$(date -u -d "${d:0:8} ${d:8:2}:${d:10:2} UTC" +%s)
    fh=$(printf %03d $(( (td - ib) / 3600 )))
    for m in $(seq -f %03g $((10#$BMI)) $((10#$BME))); do
      f=$(ls $BDYPATH/gefs.${INI_BDY_TAG:0:8}/${INI_BDY_TAG:8:2}/$BDYPREFIX/$m/g*pgrb2*.f$fh 2>/dev/null | head -1)
      if [ -n "$f" ] && [ -r "$f" ]; then nfiles=$((nfiles+1)); else
        echo "        missing $BDYPATH/gefs.${INI_BDY_TAG:0:8}/${INI_BDY_TAG:8:2}/$BDYPREFIX/$m/g*pgrb2*.f$fh"; nmiss=$((nmiss+1)); fi
    done
  done
  check "GEFS files for the forecast ($nfiles found)" test $nmiss -eq 0
  ;;
links)
  nmiss=0
  for d in $BDY_DATES; do for m in $MEMS; do
    f=$(ls $BDYARG/*/$m/met_em.d01.$(wpsname $d).nc 2>/dev/null | head -1)
    [ -n "$f" ] || { echo "        missing met_em.d01.$(wpsname $d).nc for member $m under $BDYARG"; nmiss=$((nmiss+1)); }
  done; done
  check "met_em to link ($(echo $BDY_DATES | wc -w) times x $(echo $MEMS | wc -w) members) under $BDYARG" test $nmiss -eq 0
  ;;
esac

echo "== Configuration changes"
if [ $BDYMODE == gefs ]; then PERTLINE="(boundary/target members and BDY_PERT of the source)"
else PERTLINE="RUN_BDY_PERT=0  BDY_PERT=0  BDY_MEM_INI=$MEM_INI  BDY_MEM_END=$MEM_END"; fi
cat <<EOF
  conf/config.env  : BASEDIR=$NEWABS/
  conf/exp.conf    : EXPTYPE=DAFCST  INI_STEP=0  RUN_LETKF=0  RUN_FCST=1  $PERTLINE
                     FCST_INI_DATE=FCST_END_DATE='$INIT'  FCST_INI_FREQ=300
                     FCST_LEAD_TIME=$LEAD  FCST_OFREQ=300
EOF
case $BDYMODE in
  src)   echo "                     RUN_WPS=0  EXTWPSPATH=$SRC/HIST/WPS/met_em" ;;
  anal)  echo "                     RUN_WPS=1  WPS_DATA_SOURCE='WRF'  BDYPATH=$NEWABS/BDY/  EXTWPSPATH=''"
         echo "  BDY/$INI_BDY_TAG/<MEM>/wrfout_d01_<date> -> $BDYARG/<date>/anal<MEM>"
         echo "  WPS/code/* -> $SRC/WPS/code/* except wrf_to_wps.exe -> $W2W"
         echo "  WPS/geogrid/geo_em.d01.nc -> $SRC/WPS/geogrid/"
         echo "  conf/model.conf  : METLEV=$((METLEV_SRC+1))" ;;
  gefs)  echo "                     RUN_WPS=1  EXTWPSPATH=''  (WPS from GEFS and pert_met_em in the job)"
         echo "  WPS/code, PERT/code: links to $SRC (namelists excluded); WPS/Vtable*, geo_em linked" ;;
  links) echo "                     RUN_WPS=0  EXTWPSPATH=$NEWABS/BDY_MET_EM"
         echo "  BDY_MET_EM/$INI_BDY_TAG/<MEM>/met_em.d01.<date>.nc -> files under $BDYARG" ;;
esac
cat <<EOF
  conf/machine.conf: INODE=$INODE  TOTAL_TIME_LIMIT=$TLIM
  HIST/ANAL/$INIT_TAG -> $ANAL_DIR
  WRF/code/*  -> $SRC/WRF/code/*
EOF

if [ $DRY -eq 1 ]; then
  echo "== Dry run: $nfail failed checks, nothing created"
  exit $nfail
fi
[ $nfail -eq 0 ] || { echo "== $nfail failed checks, nothing created"; exit 1; }

echo "== Building $NEWABS"
mkdir -p $NEWABS
NEW=$NEWABS
#Scripts and configuration, without job logs, caches and scratch files
#(some of them are not readable by other users of the group).
( cd $SRC && tar -cf - $TAR_EXCLUDE bin lib conf NAMELIST ) | ( cd $NEW && tar -xf - )
rm -f $NEW/conf/step.conf
mkdir -p $NEW/WRF/code $NEW/HIST/ANAL $NEW/HIST/DAFCST $NEW/LOGS $NEW/PROCS
ln -s $SRC/WRF/code/* $NEW/WRF/code/
ln -s $ANAL_DIR $NEW/HIST/ANAL/$INIT_TAG
stripe () { if command -v lfs > /dev/null ; then lfs setstripe -c ${2:-8} $1 ; fi ; }
stripe $NEW/HIST/DAFCST
stripe $NEW/WRF 4          #wrfinput, wrfbdy and interpolated met_em of every member

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
setvar $NEW/conf/exp.conf    INI_STEP        0      #the source may hold a restart step
setvar $NEW/conf/exp.conf    RUN_LETKF       0
setvar $NEW/conf/exp.conf    RUN_FCST        1
if [ $BDYMODE != gefs ]; then
  setvar $NEW/conf/exp.conf  RUN_BDY_PERT    0
  setvar $NEW/conf/exp.conf  BDY_PERT        0
  setvar $NEW/conf/exp.conf  BDY_MEM_INI     $MEM_INI
  setvar $NEW/conf/exp.conf  BDY_MEM_END     $MEM_END
fi
setvar $NEW/conf/exp.conf    FCST_INI_DATE   "'$INIT'"
setvar $NEW/conf/exp.conf    FCST_END_DATE   "'$INIT'"
setvar $NEW/conf/exp.conf    FCST_INI_FREQ   300
setvar $NEW/conf/exp.conf    FCST_LEAD_TIME  $LEAD
setvar $NEW/conf/exp.conf    FCST_OFREQ      300
setvar $NEW/conf/machine.conf INODE          $INODE
setvar $NEW/conf/machine.conf TOTAL_TIME_LIMIT "\"$TLIM\""

case $BDYMODE in
src)
  setvar $NEW/conf/exp.conf RUN_WPS    0
  setvar $NEW/conf/exp.conf EXTWPSPATH "\"$SRC/HIST/WPS/met_em\""
  ;;
anal)
  setvar $NEW/conf/exp.conf RUN_WPS         1
  setvar $NEW/conf/exp.conf WPS_DATA_SOURCE "'WRF'"
  setvar $NEW/conf/exp.conf BDYPATH         "\"$NEW/BDY/\""
  setvar $NEW/conf/exp.conf EXTWPSPATH      '""'
  for m in $MEMS; do
    mkdir -p $NEW/BDY/$INI_BDY_TAG/$m
    for d in $BDY_DATES; do
      ln -s $BDYARG/$d/anal00$m $NEW/BDY/$INI_BDY_TAG/$m/wrfout_d01_$(wpsname $d)
    done
  done
  mkdir -p $NEW/WPS/code $NEW/WPS/geogrid $NEW/HIST/WPS/met_em_ori
  ln -s $SRC/WPS/code/* $NEW/WPS/code/
  rm -f $NEW/WPS/code/wrf_to_wps.exe
  cp $W2W $NEW/WPS/code/wrf_to_wps.exe      #a copy, so a later rebuild cannot change a running experiment
  setvar $NEW/conf/model.conf METLEV $((METLEV_SRC+1))
  ln -s $SRC/WPS/geogrid/geo_em.d01.nc $NEW/WPS/geogrid/
  stripe $NEW/HIST/WPS/met_em_ori 4
  stripe $NEW/WPS 4         #intermediate files of wrf_to_wps
  ;;
gefs)
  setvar $NEW/conf/exp.conf RUN_WPS    1
  setvar $NEW/conf/exp.conf EXTWPSPATH '""'
  for d in WPS/code PERT/code; do
    mkdir -p $NEW/$d
    for f in $SRC/$d/*; do
      case $(basename $f) in *namelist*) continue ;; esac   #the run scripts cp over namelists
      ln -s $(readlink -f $f) $NEW/$d/
    done
  done
  mkdir -p $NEW/WPS/geogrid $NEW/HIST/WPS $NEW/PERT
  ln -s $(readlink -f $SRC/WPS/geogrid/geo_em.d01.nc) $NEW/WPS/geogrid/
  for f in $SRC/WPS/Vtable*; do ln -s $(readlink -f $f) $NEW/WPS/; done
  stripe $NEW/HIST/WPS 4
  ;;
links)
  setvar $NEW/conf/exp.conf RUN_WPS    0
  setvar $NEW/conf/exp.conf EXTWPSPATH "\"$NEW/BDY_MET_EM\""
  for m in $MEMS; do
    mkdir -p $NEW/BDY_MET_EM/$INI_BDY_TAG/$m
    for d in $BDY_DATES; do
      f=$(ls $BDYARG/*/$m/met_em.d01.$(wpsname $d).nc | head -1)
      ln -s $(readlink -f $f) $NEW/BDY_MET_EM/$INI_BDY_TAG/$m/
    done
  done
  ;;
esac

#Cleanup helper for after the job: work directories only, never outputs.
cat > $NEW/clean_work.sh <<EOF
#!/bin/bash
# Remove the work directories of this experiment (made by make_fcst_exp.sh).
# Keeps HIST/ (forecasts, met_em_ori) and LOGS/.
cd $NEW || exit 1
rm -rf WRF/[0-9][0-9][0-9] WPS/[0-9][0-9][0-9] PERT/[0-9][0-9] PROCS/*
echo "Removed work directories under $NEW/WRF and $NEW/WPS"
EOF
chmod +x $NEW/clean_work.sh

echo "== Changed lines"
grep -h "make_fcst_exp.sh" $NEW/conf/*.conf $NEW/conf/config.env | sed 's/   # set by.*//; s/^/  /'
echo "== Done. Submit with: cd $NEW/bin && ./queue_sub.sh main_FCST.sh"
echo "   After the job: $NEW/clean_work.sh"
