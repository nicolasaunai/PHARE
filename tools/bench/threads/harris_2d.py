#!/usr/bin/env python3
"""
2D Harris sheet used to compare builds and thread counts: same input, timed advance loop.

    PHARE_THREADS=24 mpirun -n 2 --map-by socket:PE=24 python3 harris_2d.py

Settings come from the environment so one script serves every run of a scan:
  BENCH_CELLS      level 0 cells, "400x200"
  BENCH_PPC        particles per cell, 100
  BENCH_STEPS      coarse time steps, 20
  BENCH_LEVELS     1 (no refinement) or 2 (one refined box over both sheets), 1
  BENCH_REFINE     "boxes": the refined box is given and never moves (default);
                   "tagging": cells are tagged and the hierarchy regridded as the run goes
  BENCH_PATCH      largest_patch_size, 25: this sets how many patches a thread can pick from
  BENCH_DIAG_DIR   where to write E, B, moments and particles at the last step, none if unset
  BENCH_DIAG_EVERY with BENCH_DIAG_DIR: E, B and moments every that many steps, no particles
  BENCH_AFFINITY   if set, every rank prints the cpus its threads may run on, after the run
  BENCH_TAG        free text copied in the BENCH line, to tell the runs of a scan apart

Prints one line "BENCH key=value ..." on rank 0 with the wall times of the initialization and of
the advance loop.
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
refine = os.environ.get("BENCH_REFINE", "boxes")
largest_patch_size = int(os.environ.get("BENCH_PATCH", 25))
diag_dir = os.environ.get("BENCH_DIAG_DIR")
diag_every = int(os.environ.get("BENCH_DIAG_EVERY", 0))

dl = (0.4, 0.4)
time_step = 0.005
final_time = steps * time_step


def config():
    L = 0.5

    refinement = {}
    if levels == 2 and refine == "tagging":
        refinement = {"refinement": "tagging", "max_nbr_levels": 2, "nesting_buffer": 1}
    elif levels == 2:
        # one box over the middle half in y, where both current sheets are, clear of the
        # periodic x boundaries
        lower = (cells[0] // 8, cells[1] // 4)
        upper = (7 * cells[0] // 8 - 1, 3 * cells[1] // 4 - 1)
        refinement = {"refinement_boxes": {"L0": {"B0": [lower, upper]}}}

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
        times = [final_time]
        if diag_every:
            times = [step * time_step for step in range(diag_every, steps + 1, diag_every)]
        for quantity in ["E", "B"]:
            ph.ElectromagDiagnostics(quantity=quantity, write_timestamps=times)
        for quantity in ["charge_density", "bulkVelocity"]:
            ph.FluidDiagnostics(quantity=quantity, write_timestamps=times)
        if not diag_every:
            ph.ParticleDiagnostics(
                quantity="domain", write_timestamps=times, population_name="protons"
            )
    return sim


def affinity():
    """The cpus each thread of this process may run on: "pinned" threads have a single one."""
    allowed = []
    for task in os.listdir("/proc/self/task"):
        try:
            with open(f"/proc/self/task/{task}/status") as status:
                for line in status:
                    if line.startswith("Cpus_allowed_list"):
                        allowed.append(line.split(":")[1].strip())
        except OSError:  # the thread ended meanwhile
            pass
    pinned = sorted(int(cpus) for cpus in allowed if cpus.isdigit())
    free = sorted(set(cpus for cpus in allowed if not cpus.isdigit()))
    return (
        f"BENCH_AFFINITY rank={cpp.mpi_rank()} threads={len(allowed)}"
        f" pinned={len(pinned)} pinned_cpus={','.join(map(str, pinned)) or 'none'}"
        f" others_may_run_on={'|'.join(free) or 'none'}"
    )


def main():
    simulator = Simulator(config(), print_one_line=False)
    start = time.perf_counter()
    simulator.initialize()
    initialize = time.perf_counter() - start

    start = time.perf_counter()
    while simulator.currentTime() < final_time - time_step / 2:
        simulator.advance()
    advance = time.perf_counter() - start

    # one write per line: the lines of several ranks must not get mixed
    if os.environ.get("BENCH_AFFINITY"):
        print(affinity() + "\n", end="", flush=True)

    if cpp.mpi_rank() == 0:
        line = " ".join(
            [
                "BENCH",
                os.environ.get("BENCH_TAG", ""),
                f"ranks={cpp.mpi_size()}",
                f"threads={os.environ.get('PHARE_THREADS', 1)}",
                f"cells={cells[0]}x{cells[1]}",
                f"ppc={ppc}",
                f"levels={levels}",
                f"refine={refine}",
                f"patch={largest_patch_size}",
                f"steps={steps}",
                f"init_s={initialize:.3f}",
                f"advance_s={advance:.3f}",
            ]
        )
        print(line + "\n", end="", flush=True)
    simulator.reset()


if __name__ == "__main__":
    main()
