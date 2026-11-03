"""Diagnose differences between two LETKF output directories.

Called by run_regression.sh only when cmp finds differing files. For each
NetCDF file it lists the variables that differ, with the number of
differing values and the largest absolute difference.

Usage: python3 compare_outputs.py DIR_A DIR_B FILE [FILE ...]
"""
import sys
from pathlib import Path

import numpy as np

try:
    from netCDF4 import Dataset
except ImportError:
    sys.exit("netCDF4 not available: install it to see which variables differ")


def compare_file(path_a, path_b):
    with Dataset(path_a) as a, Dataset(path_b) as b:
        names = sorted(set(a.variables) | set(b.variables))
        n_diff = 0
        for name in names:
            if name not in a.variables or name not in b.variables:
                print(f"    {name}: present in only one file")
                n_diff += 1
                continue
            va = np.ma.filled(a.variables[name][:], np.nan)
            vb = np.ma.filled(b.variables[name][:], np.nan)
            if va.shape != vb.shape:
                print(f"    {name}: shapes {va.shape} vs {vb.shape}")
                n_diff += 1
                continue
            if va.dtype.kind in "fc":
                same = (va == vb) | (np.isnan(va) & np.isnan(vb))
            else:
                same = va == vb
            if not np.all(same):
                bad = ~same
                if va.dtype.kind in "fciu":
                    dmax = np.nanmax(np.abs(va[bad].astype(float) - vb[bad].astype(float)))
                    print(f"    {name}: {bad.sum()} of {bad.size} values differ, max|diff| = {dmax:.3e}")
                else:
                    print(f"    {name}: {bad.sum()} of {bad.size} values differ")
                n_diff += 1
        if n_diff == 0:
            print("    all variables identical (difference in header or attributes only)")


def main():
    dir_a, dir_b = Path(sys.argv[1]), Path(sys.argv[2])
    for name in sys.argv[3:]:
        fa, fb = dir_a / name, dir_b / name
        if fa.read_bytes() == fb.read_bytes():
            continue
        print(f"  {name}:")
        compare_file(fa, fb)


if __name__ == "__main__":
    main()
