#!/bin/bash
#=======================================================================
# Build the 4 km, 5-min wrfutil DA-cycle experiment on Hydra (CIMA) from
# the template osse/conf/PREVENIR_4KM_5MIN_HYDRA_2024031900, with the
# wrfutil scripts of a wrfutil_hydra clone (bin/CreateExp.sh). It is the
# Hydra counterpart of osse/make_da_exp.sh, which copied the director's
# experiment on Fugaku; on Hydra there is no experiment to copy.
#
# The template reproduces the director's Fugaku run
# PREVENIR_LOWRESOLUTION_DA_FUGAKU_2024031900_5MIN_WINIT:
#   exp.conf, letkf.conf, obs.conf  from the 5MIN_WINIT_SINGLECONF copy of
#                                   its conf (same DA settings; the
#                                   single-physics option is in model.conf)
#   model.conf, machine.conf, config.env  from wrfutil's
#                                   PREVENIR_LOWRESOLUTION_DA_HYDRA_2024031900
#                                   (multi-physics, Hydra paths), with
#                                   48-core nodes and the process layouts
#                                   of the Fugaku run
# The wrfutil clone is not modified: the template and the scripts are
# copied to a temporary folder and CreateExp.sh runs there.
#
# Usage:
#   osse/make_exp_hydra.sh [options] NEW_EXP ["DA_INI_DATE" "DA_END_DATE" [SPIN_UP_LENGTH]]
#
#   NEW_EXP  new experiment directory (must not exist), in the scratch
#            area (~/salidas = /nfsmounts/storage/scratch/$USER)
#   dates    default: the template's, 2024-03-19 06:00 -> 18:00 UTC with a
#            6 h spin-up (option B: 73 steps, analyses every 5 min 12-18 UTC)
#
# Options:
#   --dry-run       check every input and print the settings; create nothing
#   --smoke         4 members (BDY 001-004 -> target 001-004), NITER=1, 4
#                   nodes, 06:00:00, last analysis 10 min after the spin-up
#                   (3 steps) unless dates are given. Also the benchmark: the
#                   spin-up step times one 6 h WRF forecast per node.
#   --inode N       nodes of the job [12; 4 with --smoke]. WRF runs one member
#                   per node, the LETKF one process per node with 48 threads.
#   --time HH:MM:SS job time limit [96:00:00; 06:00:00 with --smoke]
#   --obspath DIR   observations [~/salidas/OBS5MIN]
#   --letkf TAR     LETKF executables [~/salidas/WRF_TEnKF_reference/dist/letkf_INTEL_HYDRA.tar,
#                   the original code built with make.sh INTEL_HYDRA]
#   --wrfutil DIR   wrfutil_hydra clone [~/salidas/wrfutil_hydra]
#
# bin/run_LETKF.sh is the WRF_TEnKF version: wrfutil's, passing the tempering
# settings of letkf.conf to the namelist (none set = the original LETKF).
# bin/run_ensfcst.sh is wrfutil's with the #WRF_TEnKF blocks (wrfinput_d01 renamed
# and the analysis copied on the member's node, with bin/wrfinput_from_anal.sh).
#
# After building: cd NEW_EXP/bin && ./queue_sub.sh main_ASIM.sh
# To restart after a time limit, set INI_STEP in NEW_EXP/conf/exp.conf to
# the first step that did not finish and submit again.
#=======================================================================
set -e
ROOT=$(cd "$(dirname "$0")/.." && pwd)
TEMPLATE=PREVENIR_4KM_5MIN_HYDRA_2024031900
TDIR=$ROOT/osse/conf/$TEMPLATE
DRY=0 ; SMOKE=0 ; TLIM="" ; INODE=""
OBSPATH=$HOME/salidas/OBS5MIN
LETKFTAR=$HOME/salidas/WRF_TEnKF_reference/dist/letkf_INTEL_HYDRA.tar
WRFUTIL=$HOME/salidas/wrfutil_hydra
while [[ "$1" == --* ]]; do
  case $1 in
    --dry-run) DRY=1 ; shift ;;
    --smoke)   SMOKE=1 ; shift ;;
    --inode)   INODE=$2 ; shift 2 ;;
    --time)    TLIM=$2 ; shift 2 ;;
    --obspath) OBSPATH=$2 ; shift 2 ;;
    --letkf)   LETKFTAR=$2 ; shift 2 ;;
    --wrfutil) WRFUTIL=$2 ; shift 2 ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
done
if [ $# -lt 1 ]; then sed -n '3,/^#====/p' "$0"; exit 1; fi
NEW=$(readlink -m "$1")
OBSPATH=$(readlink -m "$OBSPATH")
LETKFTAR=$(readlink -m "$LETKFTAR")
WRFUTIL=$(readlink -m "$WRFUTIL")
[ -d $WRFUTIL/wrfutilV3.1 ] && WRFUTIL=$WRFUTIL/wrfutilV3.1
nfail=0

check () {   # $1 = description, $2... = test
  local d=$1; shift
  if "$@" ; then echo "  ok    $d"; else echo "  FAIL  $d"; nfail=$((nfail+1)); fi
}
confval () { # $1 = file, $2 = variable: first value assigned (quotes removed)
  grep -E "^ *(export +)?$2=" $1 | head -1 | sed -E "s/^ *(export +)?$2=//; s/[[:space:]]+#.*//; s/[[:space:]]+$//; s/^['\"]//; s/['\"]$//"
}
ensfcst_script_ok () {  # $1 = this repository's run_ensfcst.sh, $2 = the source's: equal, or equal once
                      # the repository's #WRF_TEnKF{ ... #WRF_TEnKF} blocks are replaced by their #orig: lines
  grep -q "wrfinput_from_anal.sh" "$1" || return 1
  cmp -s "$1" "$2" && return 0
  diff -q <(awk '/^#WRF_TEnKF\{/{inb=1; next} /^#WRF_TEnKF\}/{inb=0; next} inb && /^#orig: /{print substr($0, 8); next} !inb' "$1") \
          "$2" > /dev/null 2>&1
}
letkf_script_ok () {  # $1 = this repository's run_LETKF.sh, $2 = wrfutil's: equal once the tempering lines are removed
  diff -q <(grep -v -e '__NTEMP__' -e '__ALPHA_S__' -e '__TEMP_RELAX_MODE__' -e 'Tempering (WRF_TEnKF)' "$1") "$2" > /dev/null 2>&1
}
setvar () {  # $1 = file, $2 = variable, $3 = value (written verbatim)
  local f=$1 v=$2 val=$3
  if grep -qE "^ *(export +)?$v=" $f ; then
    sed -i -E "s|^ *(export +)?$v=.*|export $v=$val   # set by WRF_TEnKF osse/make_exp_hydra.sh|" $f
  else
    echo "export $v=$val   # set by WRF_TEnKF osse/make_exp_hydra.sh" >> $f
  fi
}

#Settings: template values unless changed here
DA_INI=${2:-$(confval $TDIR/exp.conf DA_INI_DATE)}
SPIN=${4:-$(confval $TDIR/exp.conf SPIN_UP_LENGTH)}
AFREQ=$(confval $TDIR/exp.conf ANALYSIS_FREQ)
AWIN=$(confval $TDIR/exp.conf ANALYSIS_WIN_END)
BDY_FREQ=$(confval $TDIR/exp.conf BDY_FREQ)
BDY_INI_FREQ=$(confval $TDIR/exp.conf BDY_INI_FREQ)
BDYPATH=$(confval $TDIR/exp.conf BDYPATH)
BDYPREFIX=$(confval $TDIR/exp.conf BDYPREFIX)
BDYVTABLE=$(confval $TDIR/exp.conf BDYVTABLE)
ICORE=$(confval $TDIR/machine.conf ICORE)
PERTSKIP=$(confval $TDIR/machine.conf PERTSKIP)
PERTPROC=$(confval $TDIR/machine.conf PERTPROC)
T_INI=$(date -u -d "$DA_INI UTC" +%s)
if [ $SMOKE -eq 1 ]; then
  BMI=001 ; BME=004 ; TMI=001 ; TME=004 ; INODE=${INODE:-4} ; TLIM=${TLIM:-06:00:00}
  DA_END=${3:-$(date -u -d @$((T_INI + SPIN + 2*AFREQ)) "+%Y-%m-%d %H:%M:%S")}
else
  BMI=$(confval $TDIR/exp.conf BDY_MEM_INI) ; BME=$(confval $TDIR/exp.conf BDY_MEM_END)
  TMI=$(confval $TDIR/exp.conf MEM_INI_TAR) ; TME=$(confval $TDIR/exp.conf MEM_END_TAR)
  INODE=${INODE:-$(confval $TDIR/machine.conf INODE)} ; TLIM=${TLIM:-$(confval $TDIR/machine.conf TOTAL_TIME_LIMIT)}
  DA_END=${3:-$(confval $TDIR/exp.conf DA_END_DATE)}
fi
#pert_met_em runs on PERTPROC cores (whole nodes when above ICORE): no more than the job has
[ $PERTPROC -gt $((INODE * ICORE)) ] && PERTPROC=$((INODE * ICORE / PERTSKIP * PERTSKIP))
T_END=$(date -u -d "$DA_END UTC" +%s)
NSTEPS=$(( (T_END - T_INI - SPIN) / AFREQ + 1 ))
tag () { date -u -d @$1 +%Y%m%d%H%M%S; }

echo "== Settings"
echo "  info  template  $TDIR"
echo "  info  wrfutil   $WRFUTIL ($(git -C $WRFUTIL describe --always --dirty 2>/dev/null))"
echo "  info  spin-up $(tag $T_INI) -> $(tag $((T_INI+SPIN))) ($SPIN s), $NSTEPS steps (0..$((NSTEPS-1))),"
echo "        last analysis $(tag $((T_INI+SPIN+(NSTEPS-1)*AFREQ)))"
NTAR=$(( 10#$TME - 10#$TMI + 1 ))
echo "  info  members: boundary $BMI-$BME, target $TMI-$TME; INODE=$INODE (x$ICORE cores), limit $TLIM"
echo "  info  WRF one member per node: $INODE at a time, $(( (NTAR + INODE - 1) / INODE )) wave(s) per guess;"
echo "        LETKF $INODE processes x $ICORE threads; pert_met_em $PERTPROC cores"

echo "== Checks"
check "SPIN_UP_LENGTH is a multiple of ANALYSIS_FREQ ($SPIN, $AFREQ)" test $((SPIN % AFREQ)) -eq 0
check "at least one DA step after the spin-up" test $NSTEPS -ge 2
check "pert_met_em cores ($PERTPROC) a multiple of PERTSKIP ($PERTSKIP)" test $((PERTPROC % PERTSKIP)) -eq 0 -a $PERTPROC -gt 0
for f in bin/CreateExp.sh bin/main_ASIM.sh bin/queue_sub.sh bin/run_LETKF.sh lib/encolarPBSSSH.sh \
         namelists/SMN/letkf.namelist namelists/SMN/namelist.input.006 vtables/$BDYVTABLE; do
  check "wrfutil $f" test -r $WRFUTIL/$f
done
check "wrfutil bin/run_ensfcst.sh = the original of WRF_TEnKF's (its #orig: lines)" ensfcst_script_ok $ROOT/bin/run_ensfcst.sh $WRFUTIL/bin/run_ensfcst.sh
check "wrfutil bin/run_LETKF.sh = WRF_TEnKF's without the tempering lines" letkf_script_ok $ROOT/bin/run_LETKF.sh $WRFUTIL/bin/run_LETKF.sh
for f in $TDIR/{config.env,exp.conf,letkf.conf,machine.conf,model.conf,obs.conf}; do
  check "template $(basename $f)" test -r $f
done
#Executables named in config.env (read after the placeholders, which are not paths)
cfg () { bash -c "source $TDIR/config.env > /dev/null 2>&1; echo \${$1}"; }
WPSPATH=$(cfg WPSPATH) ; WRFPATH=$(cfg WRFPATH) ; WRFDAPATH=$(cfg WRFDAPATH) ; GEOG=$(cfg GEOG)
WRFTOWPSPATH=$(cfg WRFTOWPSPATH) ; PERTMETEMPATH=$(cfg PERTMETEMPATH)
for f in $WPSPATH/geogrid.exe $WPSPATH/ungrib.exe $WPSPATH/metgrid.exe $WRFPATH/run/real.exe \
         $WRFPATH/run/wrf.exe $WRFDAPATH/var/da/da_update_bc.exe $WRFTOWPSPATH $PERTMETEMPATH $LETKFTAR; do
  check "readable $f" test -r $f
done
check "LETKF tar holds letkf.exe and update_wrf_time.exe" bash -c "tar tf $LETKFTAR 2>/dev/null | grep -qx letkf.exe && tar tf $LETKFTAR | grep -qx update_wrf_time.exe"
check "GEOG data $GEOG" test -d $GEOG
check "/opt/load-libs.sh (machine.conf, MACHINE=HYDRA)" test -r /opt/load-libs.sh
check "new experiment does not exist yet" test ! -e $NEW
case $(readlink -m $(dirname $NEW))/ in
  /nfsmounts/storage/scratch/*) check "new experiment in the scratch area" true ;;
  *) if [ -n "$ALLOW_ANY_DIR" ]; then check "new experiment outside the scratch area (ALLOW_ANY_DIR, tests only)" true
     else check "new experiment in the scratch area (got $NEW)" false ; fi ;;
esac

#GEFS files needed by every step (lib/timeutils.sh with WPS_CYCLE=1), as in make_da_exp.sh
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

#Observations at every analysis time
RADARS=$(bash -c "source $TDIR/obs.conf > /dev/null 2>&1; echo \${RADARC_LIST[@]}")
OBSTYPES=$(bash -c "source $TDIR/obs.conf > /dev/null 2>&1; echo \${OBS_LIST[@]}")
nwith=0
for (( s=1; s<NSTEPS; s++ )); do
  ct=$(tag $((T_INI+SPIN+s*AFREQ)))
  ls $OBSPATH/RADARC/WRF_RADARC_*_$ct.dat > /dev/null 2>&1 && nwith=$((nwith+1))
done
check "observation folder $OBSPATH" test -d $OBSPATH/RADARC
echo "  info  radar files (radars: $RADARS): $nwith of $((NSTEPS-1)) analysis times have at least one"
echo "  info  observation types in obs.conf: $OBSTYPES; present: $(cd $OBSPATH 2>/dev/null && ls -d */ | tr -d / | tr '\n' ' ')"

echo "== Configuration of the new experiment (template values otherwise)"
cat <<EOF
  conf/config.env  : BASEDIR=$NEW/  LETKFPATH=$LETKFTAR
  conf/obs.conf    : OBSPATH="$OBSPATH/"
  conf/exp.conf    : DA_INI_DATE='$DA_INI'  DA_END_DATE='$DA_END'  SPIN_UP_LENGTH=$SPIN
                     BDY_MEM_INI=$BMI  BDY_MEM_END=$BME  MEM_INI_TAR=$TMI  MEM_END_TAR=$TME$( [ $SMOKE -eq 1 ] && echo "  NITER=1")
  conf/machine.conf: INODE=$INODE  TOTAL_TIME_LIMIT="$TLIM"  PERTPROC=$PERTPROC
EOF

if [ $DRY -eq 1 ]; then
  echo "== Dry run: $nfail failed checks, nothing created"
  exit $nfail
fi
[ $nfail -eq 0 ] || { echo "== $nfail failed checks, nothing created"; exit 1; }

echo "== Building $NEW"
TMP=$(mktemp -d)
trap "rm -rf $TMP" EXIT
cp -r $WRFUTIL/bin $WRFUTIL/lib $WRFUTIL/namelists $WRFUTIL/vtables $TMP/
mkdir -p $TMP/conf
cp -r $TDIR $TMP/conf/$TEMPLATE
C=$TMP/conf/$TEMPLATE
sed -i -e "s|__LETKFPATH__|$LETKFTAR|" $C/config.env
sed -i -e "s|__OBSPATH__|$OBSPATH/|" $C/obs.conf
setvar $C/exp.conf     DA_INI_DATE     "'$DA_INI'"
setvar $C/exp.conf     DA_END_DATE     "'$DA_END'"
setvar $C/exp.conf     SPIN_UP_LENGTH  $SPIN
setvar $C/exp.conf     BDY_MEM_INI     $BMI
setvar $C/exp.conf     BDY_MEM_END     $BME
setvar $C/exp.conf     MEM_INI_TAR     $TMI
setvar $C/exp.conf     MEM_END_TAR     $TME
#pert_met_em computes the members in NITER batches and runs past the last
#member when NITER exceeds the number of members (see make_da_exp.sh)
[ $SMOKE -eq 1 ] && setvar $C/exp.conf NITER 1
setvar $C/machine.conf INODE            $INODE
setvar $C/machine.conf TOTAL_TIME_LIMIT "\"$TLIM\""
setvar $C/machine.conf PERTPROC         $PERTPROC
( cd $TMP/bin && bash CreateExp.sh $TEMPLATE $NEW ) > $TMP/create.log 2>&1 || { cat $TMP/create.log; exit 1; }
mkdir -p $NEW/HIST/ANAL $NEW/HIST/GUES $NEW/HIST/WPS
cp $TMP/create.log $NEW/LOGS/create_exp.log
#run_LETKF.sh of WRF_TEnKF instead of wrfutil's: wrfutil's, passing the
#tempering settings (none set: the original LETKF).
cp $ROOT/bin/run_LETKF.sh $NEW/bin/run_LETKF.sh
#run_ensfcst.sh of WRF_TEnKF: renames and copies on the member's node (stale NFS names)
cp -p $ROOT/bin/run_ensfcst.sh $ROOT/bin/wrfinput_from_anal.sh $NEW/bin/

echo "== Checks of the new experiment"
nfail=0
check "BASEDIR in conf/config.env" grep -q "^export BASEDIR=$NEW " $NEW/conf/config.env
for f in WPS/wps.tar WPS/wrf_to_wps.tar WPS/$BDYVTABLE WRF/wrf.tar WRF/wrfda.tar WRF/pert_met_em.tar \
         PERT/pert_met_em.tar LETKF/letkf.tar NAMELIST/letkf.namelist bin/main_ASIM.sh lib/encolarPBSSSH.sh; do
  check "$f" test -s $NEW/$f
done
check "no placeholders left in conf/" bash -c "! grep -l '__[A-Z]*__' $NEW/conf/*"
check "bin/run_LETKF.sh = WRF_TEnKF's" cmp -s $ROOT/bin/run_LETKF.sh $NEW/bin/run_LETKF.sh
check "bin/run_ensfcst.sh renames wrfinput_d01 on the member's node (WRF_TEnKF)" grep -q "wrfinput_from_anal.sh" $NEW/bin/run_ensfcst.sh
check "bin/wrfinput_from_anal.sh executable" test -x $NEW/bin/wrfinput_from_anal.sh
echo "== Changed lines"
grep -h "make_exp_hydra.sh" $NEW/conf/*.conf | sed 's/   # set by.*//; s/^/  /'
[ $nfail -eq 0 ] || { echo "== $nfail failed checks in $NEW"; exit 1; }
echo "== Done. Submit with: cd $NEW/bin && ./queue_sub.sh main_ASIM.sh"
