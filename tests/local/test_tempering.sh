#!/bin/bash
#=======================================================================
# Local tests of the tempering loop on synthetic wrfout files
# (needs mpif90, NetCDF-Fortran, python3 with numpy and netCDF4).
#
#   A. ntemp = 1 is bit for bit the code before the tempering loop
#      (git revision BASE, default 0904e58, itself checked bit for bit
#      against the reference build on Fugaku), with radar observations and
#      with inflation + RTPP and inflation + RTPS, which exercise the
#      changed code paths.
#   B. Linear check: temperature observations, a localization much wider
#      than the domain, no inflation, no relaxation. In this limit Nt
#      tempered LETKF steps give the single-step analysis exactly, since
#      (I + a1 A)^(-1/2) (I + a2 A (I + a1 A)^(-1))^(-1/2) = (I + A)^(-1/2)
#      when a1 + a2 = 1. ntemp = 2 and 4 must match ntemp = 1 for T.
#   C. Reflectivity with ntemp = 3: runs, recomputes H(x) at every step,
#      reduces the QRAIN error; RTPP in both modes runs.
#
# Usage: [KEEP=1] tests/local/test_tempering.sh [BASE]   (KEEP=1 keeps the work directory)
#=======================================================================
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
BASE=${1:-0904e58}
W=$(mktemp -d)
trap "cd /; git -C $ROOT worktree remove --force $W/base > /dev/null 2>&1; if [ -n \"$KEEP\" ]; then echo work directory: $W; else rm -rf $W; fi" EXIT
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ulimit -s unlimited
MPI="mpirun --oversubscribe"

echo "== builds"
( cd $ROOT && ./make.sh GNU_LOCAL > $W/build_new.log 2>&1 )
NEW=$ROOT/build/GNU_LOCAL/letkf.exe
git -C $ROOT worktree add -q $W/base $BASE
( cd $W/base && ./make.sh GNU_LOCAL > $W/build_base.log 2>&1 )
OLD=$W/base/build/GNU_LOCAL/letkf.exe
echo "  new: $(git -C $ROOT describe --always --dirty), base: $BASE"

python3 -I $ROOT/tests/local/fake_wrfout.py $W/wrf 5
mkdir $W/syn && cd $W/syn
for m in 1 2 3 4; do ln -s ../wrf/mem0000$m gs010000$m; done
ln -s ../wrf/mem00005 truth
printf '&synth_obs\n  nens = 4\n  seed = 60\n/\n' > synth_obs.namelist
$MPI -n 2 $ROOT/build/GNU_LOCAL/synth_obs.exe > log 2>&1

#Temperature observations (id 3073) at 900, 800, 700 hPa, from the truth.
python3 -I - $W/wrf/mem00005 $W/tobs.dat <<'EOF'
import sys, numpy as np
from netCDF4 import Dataset
nc = Dataset(sys.argv[1])
p = nc['P'][0] + nc['PB'][0]
tk = (nc['T'][0] + 300.0) * (p / 1.0e5) ** (287.0 / 1005.7)
lat, lon = nc['XLAT'][0], nc['XLONG'][0]
recs = []
for (i, j) in [(13, 10), (15, 12), (11, 8), (20, 15), (8, 6)]:
    for plev in (900.0, 800.0, 700.0):
        k = np.argmin(np.abs(p[:, j, i] - plev * 100))
        recs.append([3073, lon[j, i], lat[j, i], plev, tk[k, j, i] + 0.3, 1.0, 1])
with open(sys.argv[2], 'wb') as f:
    for r in recs:
        b = np.array(r, '>f4').tobytes()
        f.write(np.int32(len(b)).astype('>i4').tobytes() + b + np.int32(len(b)).astype('>i4').tobytes())
print(f"  {len(recs)} temperature observations")
EOF

# run DIR EXE NAMELIST_TEXT OBS(radar|temp)
run () {
  local d=$1 exe=$2 nml=$3 obs=$4
  mkdir -p $d && cd $d
  for m in 1 2 3 4; do cp $W/wrf/mem0000$m gs010000$m; done
  cp $W/wrf/mem00004 guesemean; cp $W/wrf/mem00004 analemean
  if [ $obs == radar ]; then ln -sf $W/syn/rad0101.dat rad0101.dat; else ln -sf $W/tobs.dat obs01.dat; fi
  printf "$nml" > letkf.namelist
  $MPI -n 2 $exe > log 2>&1 || { echo "  letkf failed in $d"; tail -20 NOUT-00000; exit 1; }
  grep -q "LETKF ends SUCCESFULLY" NOUT-00000
}
cap4 () { # $1 = extra inflation lines, $2 = tempering section or ""
  sed -e 's/__NBV__/4/' -e 's/__NRADAR__/1/' $ROOT/namelists/letkf_osse_cap4.namelist | \
    awk -v infl="$1" -v temp="$2" '
      /^&tempering/ {skip=1}
      skip && /^\// {skip=0; if (temp != "") print temp; next}
      skip {next}
      /^&inflation/ {print; print infl; next}
      /cov_infl_mul|relax_alpha/ && infl != "" {next}
      {print}' | sed 's/%/%%/g'
}

echo "== A. ntemp = 1 against $BASE (bit for bit)"
for case in plain rtpp rtps; do
  case $case in
    plain) INFL="" ;;
    rtpp)  INFL=$'  cov_infl_mul = 1.2d0\n  relax_alpha = 0.6d0\n  relax_alpha_spread = 0.0d0' ;;
    rtps)  INFL=$'  cov_infl_mul = 1.2d0\n  relax_alpha = 0.0d0\n  relax_alpha_spread = 0.8d0' ;;
  esac
  NML=$(cap4 "$INFL" "")
  run $W/A_$case/old $OLD "$NML" radar
  run $W/A_$case/new $NEW "$NML" radar
  for f in gs0100001 gs0100002 gs0100003 gs0100004 analemean guesemean; do
    cmp -s $W/A_$case/old/$f $W/A_$case/new/$f || { echo "  FAIL $case: $f differs"; exit 1; }
  done
  echo "  ok    $case: analyses and means identical"
done

echo "== B. linear limit: temperature obs, global localization"
LIN=$'&general\n nslots=1\n nbslot=1\n nbv=4\n/\n&localization\n sigma_obs=1.0d8\n sigma_obsv=1.0d3\n lev_update_q=90\n/\n&observations\n gross_error=100.0d0\n/\n&inflation\n cov_infl_mul=1.0d0\n relax_alpha=0.0d0\n relax_alpha_spread=0.0d0\n/\n&radar_da\n nradar=0\n/\n'
run $W/B1 $NEW "$LIN&tempering\n ntemp=1\n/\n" temp
for nt in 2 4; do
  run $W/B$nt $NEW "$LIN&tempering\n ntemp=$nt\n alpha_s=1.0d0\n/\n" temp
  grep -c "TEMPERING STEP" $W/B$nt/NOUT-00000 | grep -qx $((nt-1)) || { echo "  FAIL: expected $((nt-1)) recomputations of H(x)"; exit 1; }
done
python3 -I - $W <<'EOF'
import sys, numpy as np
from netCDF4 import Dataset
W = sys.argv[1]
tb = np.array([Dataset(f'{W}/wrf/mem0000{m}')['T'][0] for m in range(1, 5)])
t1 = np.array([Dataset(f'{W}/B1/gs010000{m}')['T'][0] for m in range(1, 5)])
inc = np.abs(t1 - tb).max()
assert inc > 0.05, f"increment too small to test: {inc}"
for nt in (2, 4):
    tn = np.array([Dataset(f'{W}/B{nt}/gs010000{m}')['T'][0] for m in range(1, 5)])
    d = np.abs(tn - t1).max()
    print(f"  ntemp={nt}: max |T(ntemp) - T(1)| = {d:.2e} K, max increment {inc:.2f} K, ratio {d/inc:.1e}")
    assert d / inc < 1e-3
EOF

echo "== C. reflectivity, ntemp = 3"
for mode in none 1 2; do
  case $mode in none) INFL="" ;; *) INFL=$'  cov_infl_mul = 1.1d0\n  relax_alpha = 0.5d0\n  relax_alpha_spread = 0.0d0' ;; esac
  TMP=$'&tempering\n  ntemp = 3\n  alpha_s = 1.0d0\n  temp_relax_mode = '${mode/none/1}$'\n/'
  run $W/C_$mode $NEW "$(cap4 "$INFL" "$TMP")" radar
  grep -E "TEMPERING STEP" $W/C_$mode/NOUT-00000 | sed "s/^/  [$mode] /"
done
grep -q "RTPP after the last step" $W/C_1/NOUT-00000
! grep -q "RTPP after the last step" $W/C_2/NOUT-00000
python3 -I - $W <<'EOF'
import sys, numpy as np
from netCDF4 import Dataset
W = sys.argv[1]
qt = Dataset(f'{W}/wrf/mem00005')['QRAIN'][0]
qb = np.mean([Dataset(f'{W}/wrf/mem0000{m}')['QRAIN'][0] for m in range(1, 5)], 0)
eb = np.sqrt(np.mean((qb - qt) ** 2))
for case in ('A_plain/new', 'C_none', 'C_1', 'C_2'):
    q = np.array([Dataset(f'{W}/{case}/gs010000{m}')['QRAIN'][0] for m in range(1, 5)])
    assert np.all(np.isfinite(q)) and q.min() >= 0.0
    ea = np.sqrt(np.mean((q.mean(0) - qt) ** 2))
    print(f"  {case:12s} QRAIN RMSE {eb:.2e} -> {ea:.2e}, spread {q.std(0, ddof=1).mean():.2e}")
    assert ea < eb
EOF
echo "ALL TEMPERING TESTS PASSED"
