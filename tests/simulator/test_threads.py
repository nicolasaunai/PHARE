#!/usr/bin/env python3
"""
An MPI process can run on several threads (PHARE_THREADS): its patches concurrently, and the
particles of a patch in pieces when it has many.

With several threads the particles of a patch are deposited as they are pushed, in pieces that
do not depend on the number of threads: 2 and 4 threads must write identical files, fields,
moments and particles, on both levels, bit for bit. This is checked with patches of few
particles (one piece, many patches per thread) and of many particles (several pieces).

One thread pushes, selects, then deposits: nodes sum the same particles in another order, and
the files agree with the threaded ones to rounding.
"""

import os
import unittest
from pathlib import Path

import h5py
import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator
from tests.simulator import SimulatorTest

ph.NO_GUI()

cells = (40, 20)
dl = (0.2, 0.2)
time_step = 0.005
time_step_nbr = 10
final_time = time_step * time_step_nbr
diag_dir = "phare_outputs/threads"


def config(diag_subdir, largest_patch_size, ppc):
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells,
        dl=dl,
        interp_order=1,
        largest_patch_size=largest_patch_size,
        refinement_boxes={"L0": {"B0": [(10, 5), (29, 14)]}},
        diag_options={
            "format": "phareh5",
            "options": {"dir": f"{diag_dir}/{diag_subdir}", "mode": "overwrite"},
        },
        strict=True,
    )

    Lx, Ly = cells[0] * dl[0], cells[1] * dl[1]

    def density(x, y):
        return 1.0 + 0.3 * np.cos(2 * np.pi * x / Lx) * np.cos(2 * np.pi * y / Ly)

    def bx(x, y):
        return 1.0 + 0.0 * x

    def by(x, y):
        return 0.2 * np.sin(2 * np.pi * x / Lx)

    def zero(x, y):
        return 0.0 * x

    def vx(x, y):
        return 0.3 * np.sin(2 * np.pi * y / Ly)

    def vth(x, y):
        return 0.3 + 0.0 * x

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=zero,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": vx,
            "vbulky": zero,
            "vbulkz": zero,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 1337},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.1)

    timestamps = [0.0, final_time]
    for quantity in ["E", "B"]:
        ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=timestamps)
    for quantity in ["charge_density", "bulkVelocity"]:
        ph.FluidDiagnostics(quantity=quantity, write_timestamps=timestamps)
    ph.ParticleDiagnostics(
        quantity="domain", write_timestamps=timestamps, population_name="protons"
    )
    return sim


def datasets(h5_path):
    """{dataset path: values} of a file"""
    found = {}
    with h5py.File(h5_path, "r") as h5:
        h5.visititems(
            lambda name, obj: found.__setitem__(name, obj[...])
            if isinstance(obj, h5py.Dataset)
            else None
        )
    return found


class ThreadsTest(SimulatorTest):
    def tearDown(self):
        super().tearDown()
        ph.global_vars.sim = None
        os.environ.pop("PHARE_THREADS", None)

    def run_with(self, nthreads, case, **kwargs):
        ph.global_vars.sim = None
        os.environ["PHARE_THREADS"] = str(nthreads)
        subdir = f"{case}_{nthreads}"
        Simulator(config(subdir, **kwargs)).run().reset()
        return Path(diag_dir) / subdir

    def files_of(self, one, other):
        files = sorted(p.name for p in one.glob("*.h5"))
        self.assertGreaterEqual(len(files), 5)
        self.assertEqual(files, sorted(p.name for p in other.glob("*.h5")))
        for name in files:
            a, b = datasets(one / name), datasets(other / name)
            self.assertGreater(len(a), 0, name)
            self.assertEqual(sorted(a), sorted(b), name)
            yield name, a, b

    def assert_same_files(self, one, other):
        for name, a, b in self.files_of(one, other):
            for key in a:
                # ghost values never filled are NaN in both runs
                same = np.array_equal(a[key], b[key], equal_nan=a[key].dtype.kind == "f")
                self.assertTrue(same, f"{name}:{key} differs")

    def assert_close_files(self, one, other, rtol):
        for name, a, b in self.files_of(one, other):
            for key in a:
                self.assertEqual(a[key].shape, b[key].shape, f"{name}:{key}")
                if a[key].dtype.kind != "f":
                    self.assertTrue(np.array_equal(a[key], b[key]), f"{name}:{key} differs")
                    continue
                scale = np.nanmax(np.abs(a[key])) if a[key].size else 0.0
                close = np.allclose(a[key], b[key], rtol=0, atol=rtol * scale, equal_nan=True)
                self.assertTrue(close, f"{name}:{key} differs by more than rounding")

    def test_small_patches_same_files_on_two_and_four_threads(self):
        # 2000 particles per patch: one piece, many patches per thread
        case = dict(largest_patch_size=10, ppc=20)
        one = self.run_with(1, "small", **case)
        two = self.run_with(2, "small", **case)
        four = self.run_with(4, "small", **case)

        if cpp.mpi_rank() == 0:
            self.assert_same_files(two, four)
            self.assert_close_files(one, four, rtol=1e-12)

    def test_large_patches_same_files_on_two_and_four_threads(self):
        # 40000 particles per patch and more: push and deposit in pieces
        case = dict(largest_patch_size=20, ppc=100)
        one = self.run_with(1, "large", **case)
        two = self.run_with(2, "large", **case)
        four = self.run_with(4, "large", **case)

        if cpp.mpi_rank() == 0:
            self.assert_same_files(two, four)
            self.assert_close_files(one, four, rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
