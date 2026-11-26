#!/bin/bash
#=======================================================================
# WRF_TEnKF: replace real.exe's wrfinput_d01 by the analysis, on the member's
# node. run_ensfcst.sh runs this with $MPIEXESERIAL, i.e. on the node where
# da_update_bc and wrf.exe then open wrfinput_d01, so that node renames and
# creates the names itself (Hydra's NFS name cache, see run_ensfcst.sh).
# The analysis is read with direct I/O, bypassing this node's page cache.
# Usage: wrfinput_from_anal.sh MEMBER_DIR ANALYSIS_FILE    (exit 1 on any failure;
#        nothing is renamed unless both files exist and are not empty)
#=======================================================================
set -e
D=$1 ; A=$2
[ -s "$D/wrfinput_d01" ] || { echo "wrfinput_from_anal.sh: no $D/wrfinput_d01" ; exit 1 ; }
[ -s "$A" ]              || { echo "wrfinput_from_anal.sh: no analysis $A" ; exit 1 ; }
mv "$D/wrfinput_d01" "$D/wrfinput_d01.org"
dd if="$A" of="$D/wrfinput_d01" iflag=direct bs=16M conv=fsync status=none
[ "$(stat -c %s "$D/wrfinput_d01")" -eq "$(stat -L -c %s "$A")" ] || { echo "wrfinput_from_anal.sh: size of the copy differs from $A" ; exit 1 ; }
