#!/bin/bash
#=======================================================================
# Control forecast without DA: link the guesses of a DA run under the
# names a forecast experiment expects as initial conditions.
#   osse/link_gues_as_anal.sh GUES_DIR NMEM OUT_DIR
# makes OUT_DIR/anal<NNNNN> -> GUES_DIR/gues<NNNNN> for members 1..NMEM.
# Then: osse/make_fcst_exp.sh --gefs SRC NEW "date" LEAD 001 <NMEM> OUT_DIR
#=======================================================================
set -e
[ $# -eq 3 ] || { sed -n '3,/^#====/p' "$0"; exit 1; }
GUES=$(readlink -f "$1") ; NMEM=$((10#$2)) ; OUT=$3
[ -e $OUT ] && { echo "$OUT exists"; exit 1; }
for m in $(seq -f %05g 1 $NMEM); do
  [ -r $GUES/gues$m ] || { echo "missing $GUES/gues$m"; exit 1; }
done
mkdir -p $OUT
for m in $(seq -f %05g 1 $NMEM); do ln -s $GUES/gues$m $OUT/anal$m; done
echo "Linked $NMEM guesses of $GUES as $OUT/anal*"
