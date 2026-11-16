#!/usr/bin/env python3
"""
Compares the diagnostics of two runs, dataset by dataset:

    compare.py <diag_dir> <other_diag_dir>

Prints the number of datasets, how many differ at all, and the worst difference relative to the
largest value of the dataset. "different 0" means the two runs wrote the same bits.
"""

import sys
from pathlib import Path

import h5py
import numpy as np


def datasets(path):
    found = {}

    def keep(name, obj):
        if isinstance(obj, h5py.Dataset):
            found[name] = obj[...]

    with h5py.File(path, "r") as h5:
        h5.visititems(keep)
    return found


def main():
    one, other = Path(sys.argv[1]), Path(sys.argv[2])
    files = sorted(p.name for p in one.glob("*.h5"))
    if not files or files != sorted(p.name for p in other.glob("*.h5")):
        print(f"FILES differ or missing: {files} / {sorted(p.name for p in other.glob('*.h5'))}")
        return 1

    n = different = 0
    worst = 0.0
    for name in files:
        a, b = datasets(one / name), datasets(other / name)
        if sorted(a) != sorted(b):
            print(f"DATASETS differ in {name}: {len(a)} / {len(b)}")
            different += 1
            continue
        for key in a:
            n += 1
            if a[key].shape != b[key].shape:
                print(f"SHAPE differs {name} {key}: {a[key].shape} / {b[key].shape}")
                different += 1
            elif a[key].dtype.kind != "f":
                different += not np.array_equal(a[key], b[key])
            elif not np.array_equal(a[key], b[key], equal_nan=True):
                different += 1
                scale = np.nanmax(np.abs(a[key])) or 1.0
                worst = max(worst, np.nanmax(np.abs(a[key] - b[key])) / scale)
    print(f"datasets {n} different {different} worst relative difference {worst:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
