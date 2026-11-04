#!/bin/bash
#=======================================================================
# Build a LETKF test case from a saved wrfutil experiment, with links
# only (nothing is copied from the experiment).
#
# The case reproduces what run_LETKF.sh leaves in $LETKFDIR/00/ for a
# single-slot (3D) analysis at DATE, except that the "prior" members are
# the ANALYSES saved in HIST/ANAL/DATE: wrfutil keeps no prior members.
# Physically meaningless (it assimilates into an analysis), but every
# executable gets exactly the same inputs, which is all a regression
# test needs.
#
#   gs01<mem>      -> EXPDIR/HIST/ANAL/DATE/anal<mem>, mem = 00001..NBV
#   guesemean,
#   analemean      -> anal<NBV> (templates; run_LETKF.sh uses the last member)
#   rad01<nn>.dat  -> OBSPATH/RADARC/WRF_RADARC_<radar>_DATE.dat, with nn the
#                     position of the radar in RADARC_LIST (missing radars
#                     leave a gap, as in run_LETKF.sh)
#   obs01.dat      =  conventional obs in OBS_LIST order, concatenated
#                     (the only file written here: it is small)
#   letkf.namelist =  EXPDIR/LETKF/00/letkf.namelist with nbv=NBV
#
# Usage: tests/make_case.sh EXPDIR DATE NBV CASE_DIR [OBSPATH]
#   DATE    YYYYMMDDhhmmss, e.g. 20240319180000
#   OBSPATH defaults to the one in EXPDIR/conf/obs.conf
#=======================================================================
set -e
if [ $# -lt 4 ]; then
  sed -n '3,/^#====/p' "$0"; exit 1
fi
EXPDIR=$(cd "$1" && pwd)
DATE=$2
NBV=$3
CASE_DIR=$4

ANALDIR=$EXPDIR/HIST/ANAL/$DATE
[ -d $ANALDIR ] || { echo "No $ANALDIR"; exit 1; }
[ -f $EXPDIR/LETKF/00/letkf.namelist ] || { echo "No $EXPDIR/LETKF/00/letkf.namelist"; exit 1; }

#OBSPATH, OBS_LIST and RADARC_LIST from the experiment, read without
#sourcing the whole configuration.
OBSCONF=$EXPDIR/conf/obs.conf
eval "$(grep -E '^ *(export +)?(OBSPATH|OBS_LIST|RADARC_LIST)=' $OBSCONF | sed 's/#.*//')"
OBSPATH=${5:-$OBSPATH}
[ -d "$OBSPATH" ] || { echo "OBSPATH=$OBSPATH does not exist"; exit 1; }

nslots=$(grep -i "^ *nslots *=" $EXPDIR/LETKF/00/letkf.namelist | sed 's/.*= *\([0-9]*\).*/\1/')
[ "$nslots" == "1" ] || { echo "Only single-slot experiments are supported (nslots=$nslots)"; exit 1; }

[ -e $CASE_DIR ] && { echo "$CASE_DIR already exists"; exit 1; }
mkdir -p $CASE_DIR
cd $CASE_DIR

echo "Members 1..$NBV from $ANALDIR"
for IMEM in $(seq 1 $NBV); do
  f=$ANALDIR/anal$(printf %05d $IMEM)
  [ -f $f ] || { echo "Missing $f"; exit 1; }
  ln -s $f gs01$(printf %05d $IMEM)
done
ln -s $ANALDIR/anal$(printf %05d $NBV) guesemean
ln -s $ANALDIR/anal$(printf %05d $NBV) analemean

echo "Observations at $DATE from $OBSPATH"
NRAD=1
for RADARC in "${RADARC_LIST[@]}"; do
  f=$OBSPATH/RADARC/WRF_RADARC_${RADARC}_${DATE}.dat
  if [ -e $f ]; then
    ln -s $f rad01$(printf %02d $NRAD).dat
    echo "  radar $NRAD $RADARC"
  else
    echo "  radar $NRAD $RADARC missing"
  fi
  NRAD=$(( NRAD + 1 ))
done
for OBS_TYPE in "${OBS_LIST[@]}"; do
  [ $OBS_TYPE == "RADARC" ] && continue
  f=$OBSPATH/$OBS_TYPE/WRF_${OBS_TYPE}_${DATE}.dat
  if [ -e $f ]; then
    cat $f >> obs01.dat
    echo "  conventional $OBS_TYPE"
  else
    echo "  conventional $OBS_TYPE missing"
  fi
done

sed -e "s/^\( *nbv *= *\)[0-9]*/\1$NBV/I" $EXPDIR/LETKF/00/letkf.namelist > letkf.namelist
grep -qi "^ *nbv *= *$NBV *$" letkf.namelist || { echo "Could not set nbv in letkf.namelist"; exit 1; }
grep -qi "^ *nradar *= *${#RADARC_LIST[@]}" letkf.namelist || \
  echo "WARNING: nradar in letkf.namelist differs from the ${#RADARC_LIST[@]} radars in RADARC_LIST"

echo "Case ready in $CASE_DIR"
