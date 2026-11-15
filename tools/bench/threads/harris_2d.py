#!/usr/bin/env python3
"""
2D Harris sheet used to compare builds and thread counts: same input, timed advance loop.

    PHARE_THREADS=24 mpirun -n 2 --map-by socket:PE=24 python3 harris_2d.py

Settings come from the environment so one script serves every run of a scan:
  BENCH_CELLS      level 0 cells, "400x200"
  BENCH_PPC        particles per cell, 100
  BENCH_STEPS      coarse time steps, 20
  BENCH_LEVELS     1 (no refinement) or 2 (one refined box around both sheets), 1
  BENCH_PATCH      largest_patch_size, 25: this sets how many patches a thread can pick from
  BENCH_DIAG_DIR   where to write E, B, moments and particles at the last step, none if unset

Prints one line "BENCH key=value ..." on rank 0 with the wall time of the advance loop.
"""

import os
import time

import numpy as np

import pyphare.pharein as ph
from pyphare import cpp
from pyphare.simulator.simulator import Simulator

ph.NO_GUI()

cells = tuple(int(n) for n in os.environ.get("BENCH_CELLS", "400x200").split("x"))
ppc = int(os.environ.get("BENCH_PPC", 100))
steps = int(os.environ.get("BENCH_STEPS", 20))
levels = int(os.environ.get("BENCH_LEVELS", 1))
largest_patch_size = int(os.environ.get("BENCH_PATCH", 25))
diag_dir = os.environ.get("BENCH_DIAG_DIR")

dl = (0.4, 0.4)
time_step = 0.005
final_time = steps * time_step


def config():
    L = 0.5

    refinement = {}
    if levels == 2:
        # one box over the middle half in y, where both current sheets are
        refinement = {
            "refinement_boxes": {
                "L0": {"B0": [(0, cells[1] // 4), (cells[0] - 1, 3 * cells[1] // 4 - 1)]}
            }
        }

    diag_options = {}
    if diag_dir:
        diag_options = {
            "diag_options": {
                "format": "phareh5",
                "options": {"dir": diag_dir, "mode": "overwrite"},
            }
        }

    sim = ph.Simulation(
        time_step=time_step,
        final_time=final_time,
        cells=cells,
        dl=dl,
        largest_patch_size=largest_patch_size,
        hyper_resistivity=0.002,
        resistivity=0.001,
        strict=True,
        **refinement,
        **diag_options,
    )

    Lx, Ly = sim.simulation_domain()

    def density(x, y):
        return (
            0.4
            + 1.0 / np.cosh((y - Ly * 0.3) / L) ** 2
            + 1.0 / np.cosh((y - Ly * 0.7) / L) ** 2
        )

    def S(y, y0, l):
        return 0.5 * (1.0 + np.tanh((y - y0) / l))

    def by(x, y):
        x0, y1, y2 = x - 0.5 * Lx, y - 0.3 * Ly, y - 0.7 * Ly
        return 0.2 * x0 * (np.exp(-(x0**2 + y1**2)) - np.exp(-(x0**2 + y2**2)))

    def bx(x, y):
        x0, y1, y2 = x - 0.5 * Lx, y - 0.3 * Ly, y - 0.7 * Ly
        perturbation = -0.2 * y1 * np.exp(-(x0**2 + y1**2)) + 0.2 * y2 * np.exp(
            -(x0**2 + y2**2)
        )
        return -1 + 2 * (S(y, Ly * 0.3, L) - S(y, Ly * 0.7, L)) + perturbation

    def bz(x, y):
        return 0.0 * x

    def vth(x, y):
        temperature = (0.7 - (bx(x, y) ** 2 + by(x, y) ** 2) * 0.5) / density(x, y)
        assert np.all(temperature > 0)
        return np.sqrt(temperature)

    ph.MaxwellianFluidModel(
        bx=bx,
        by=by,
        bz=bz,
        protons={
            "charge": 1,
            "density": density,
            "vbulkx": bz,
            "vbulky": bz,
            "vbulkz": bz,
            "vthx": vth,
            "vthy": vth,
            "vthz": vth,
            "nbr_part_per_cell": ppc,
            "init": {"seed": 12334},
        },
    )
    ph.ElectronModel(closure="isothermal", Te=0.0)

    if diag_dir:
        last = [final_time]
        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=last)
        for quantity in ["charge_density", "bulkVelocity"]:
            ph.FluidDiagnostics(quantity=quantity, write_timestamps=last)
        ph.ParticleDiagnostics(
            quantity="domain", write_timestamps=last, population_name="protons"
        )
    return sim


def main():
    simulator = Simulator(config(), print_one_line=False)
    simulator.initialize()

    start = time.perf_counter()
    while simulator.currentTime() < final_time - time_step / 2:
        simulator.advance()
    advance = time.perf_counter() - start

    if cpp.mpi_rank() == 0:
        print(
            "BENCH",
            f"ranks={cpp.mpi_size()}",
            f"threads={os.environ.get('PHARE_THREADS', 1)}",
            f"cells={cells[0]}x{cells[1]}",
            f"ppc={ppc}",
            f"levels={levels}",
            f"patch={largest_patch_size}",
            f"steps={steps}",
            f"advance_s={advance:.3f}",
            flush=True,
        )
    simulator.reset()


if __name__ == "__main__":
    main()
