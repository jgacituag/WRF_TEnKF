#!/bin/bash
#=======================================================================
# Local smoke test of the OSSE chain on synthetic wrfout files
# (needs ./make.sh GNU_LOCAL, python3 with numpy and netCDF4):
#   1. fake_wrfout.py writes 4 members and a truth (40x30x20 grid).
#   2. synth_obs.exe writes the synthetic reflectivity.
#   3. obs_ope.exe reads the noise-free file with radar_obs_lonlatz and
#      must recover h(x_true) at the same points.
#   4. letkf.exe assimilates the noisy file with the Chapter 4 settings
#      (namelists/letkf_osse_cap4.namelist) and must reduce the QRAIN error.
#=======================================================================
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
BIN=$ROOT/build/GNU_LOCAL
W=$(mktemp -d)
trap "rm -rf $W" EXIT
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ulimit -s unlimited
MPI="mpirun --oversubscribe"

python3 -I $ROOT/tests/local/fake_wrfout.py $W/wrf 5

echo "== synth_obs"
mkdir $W/syn && cd $W/syn
for m in 1 2 3 4; do ln -s ../wrf/mem0000$m gs010000$m; done
ln -s ../wrf/mem00005 truth
printf '&synth_obs\n  nens = 4\n  seed = 60\n/\n' > synth_obs.namelist
$MPI -n 2 $BIN/synth_obs.exe > log 2>&1
grep -q "SYNTH_OBS ends SUCCESFULLY" SYNO-00000
grep -E "Candidate|Written" SYNO-00000

echo "== obs_ope round trip"
mkdir $W/ope && cd $W/ope
ln -s ../wrf/mem00005 gs0100001
ln -s ../syn/truth_rad0101.dat rad0101.dat
printf '&general\n nslots=1\n nbslot=1\n nbv=1\n enable_checks=.false.\n/\n&radar_da\n nradar=1\n radar_obs_lonlatz=.true.\n/\n' > letkf.namelist
$MPI -n 1 $BIN/obs_ope.exe > log 2>&1
python3 -I - <<'EOF'
import numpy as np
def recs(path):
    raw = open(path, 'rb').read(); out = []; p = 0
    while p < len(raw):
        n = int(np.frombuffer(raw[p:p+4], '>i4')[0]); p += 4
        out.append(raw[p:p+n]); p += n + 4
    return out
t = np.array([np.frombuffer(r, '>f4') for r in recs('../syn/truth_rad0101.dat')[3:]])
ro = recs('obsope00001')
o = np.array([np.frombuffer(r, '>f8') for r in ro[1:]])
key = lambda a, b, c: [f"{x:.3f}_{y:.3f}_{z:.0f}" for x, y, z in zip(a, b, c)]
tz = np.where(t[:, 4] > 1.0, 10*np.log10(np.maximum(t[:, 4], 1e-30)), 0.0)
d = dict(zip(key(t[:, 1], t[:, 2], t[:, 3]), tz))
ref = np.array([d.get(x, np.nan) for x in key(o[:, 1], o[:, 2], o[:, 3])])
err = np.nanmax(np.abs(o[:, 4] - o[:, 10] - ref))
print(f"  {np.isfinite(ref).sum()} of {len(t)} points recovered, max |dH| = {err:.2e} dBZ")
assert np.isfinite(ref).sum() == len(t) == len(o) and err < 0.01
EOF

echo "== letkf with Chapter 4 settings"
mkdir $W/letkf && cd $W/letkf
for m in 1 2 3 4; do cp ../wrf/mem0000$m gs010000$m; done
cp ../wrf/mem00004 guesemean; cp ../wrf/mem00004 analemean
ln -s ../syn/rad0101.dat rad0101.dat
sed -e 's/__NBV__/4/' -e 's/__NRADAR__/1/' -e 's/__NTEMP__/1/' -e 's/__ALPHA_S__/1.0d0/' -e 's/__TEMP_RELAX_MODE__/1/' \
    $ROOT/namelists/letkf_osse_cap4.namelist > letkf.namelist
$MPI -n 2 $BIN/letkf.exe > log 2>&1
grep -q "LETKF ends SUCCESFULLY" NOUT-00000
grep -E "OBSERVATIONS TO BE ASSIMILATED|LOW RAIN|USELESS" NOUT-00000 | sed 's/^/  /'
python3 -I - <<'EOF'
import numpy as np
from netCDF4 import Dataset
qb = np.mean([Dataset(f'../wrf/mem0000{m}')['QRAIN'][0] for m in range(1, 5)], 0)
qa = np.mean([Dataset(f'gs010000{m}')['QRAIN'][0] for m in range(1, 5)], 0)
qt = Dataset('../wrf/mem00005')['QRAIN'][0]
eb, ea = np.sqrt(np.mean((qb-qt)**2)), np.sqrt(np.mean((qa-qt)**2))
print(f"  QRAIN RMSE prior {eb:.2e} -> analysis {ea:.2e}")
assert ea < eb
EOF
echo "ALL OSSE TOOL TESTS PASSED"
