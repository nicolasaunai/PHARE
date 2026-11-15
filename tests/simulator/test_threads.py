#!/usr/bin/env python3
"""
The patches of an MPI process can run on several threads (PHARE_THREADS). A patch is computed
the same way whichever thread runs it, so the same run on 1 and on 4 threads must write
identical files: fields, moments and particles, on both levels, bit for bit.
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


def config(diag_subdir):
    sim = ph.Simulation(
        time_step=time_step,
        time_step_nbr=time_step_nbr,
        cells=cells,
        dl=dl,
        interp_order=1,
        largest_patch_size=10,  # many patches per thread
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
            "nbr_part_per_cell": 20,
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

    def run_with(self, nthreads):
        ph.global_vars.sim = None
        os.environ["PHARE_THREADS"] = str(nthreads)
        Simulator(config(f"threads_{nthreads}")).run().reset()

    def test_same_files_on_one_and_four_threads(self):
        self.run_with(1)
        self.run_with(4)

        if cpp.mpi_rank() != 0:
            return

        one, four = (Path(diag_dir) / f"threads_{n}" for n in (1, 4))
        files = sorted(p.name for p in one.glob("*.h5"))
        self.assertGreaterEqual(len(files), 5)
        self.assertEqual(files, sorted(p.name for p in four.glob("*.h5")))

        for name in files:
            a, b = datasets(one / name), datasets(four / name)
            self.assertGreater(len(a), 0, name)
            self.assertEqual(sorted(a), sorted(b), name)
            for key in a:
                # ghost values never filled are NaN in both runs
                same = np.array_equal(a[key], b[key], equal_nan=a[key].dtype.kind == "f")
                self.assertTrue(same, f"{name}:{key} differs")


if __name__ == "__main__":
    unittest.main()
