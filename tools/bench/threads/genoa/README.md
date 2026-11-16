# Thread benchmarks on one Adastra GENOA node

Instructions for whoever runs them (a person or their coding agent). Everything needed is in
this directory. Nothing here has been run on Adastra yet: the scripts were only tried on a
workstation with `mpirun`, at toy size.

## What this is for

The branch `threads/patch-level` lets one MPI process of PHARE use several CPU threads
(`PHARE_THREADS`). On a 48-core machine, a few processes with many threads each take about 40 %
less time than master with one process per core. The question is whether this holds on a GENOA
node (2 x 96 cores), and which split of the 192 cores into processes x threads is best there.

All runs use one node and the same input: a 2D Harris sheet (`../harris_2d.py`), 1600 x 800
cells, 100 particles per cell (128 million particles), patches of 50 cells.

## Ground rules

- One node at a time, exclusive. The whole sequence is about 2 hours of one node; do not use
  more than 6 node-hours in total without asking the person you work for.
- Do not modify PHARE's sources. If a script of this directory needs a change to work on
  Adastra (launch line, modules), change it, and say exactly what you changed in the report.
- Do not push, publish or share anything without the agreement of the person you work for.
- Report what happened, including failures, with the log. Do not rerun until numbers look
  good, and do not drop runs that look odd.
- The SLURM account, the scratch directory and the modules are those of the person you work
  for: ask them, do not guess.

## Step 1: two builds

Two source trees, built the same way (same modules, same compiler, same options):

| name   | what                                                    |
|--------|---------------------------------------------------------|
| branch | branch `threads/patch-level` (it contains this directory) |
| master | commit `dc62fdd4e`, the commit the branch starts from   |

    git clone --recursive -b threads/patch-level https://github.com/nicolasaunai/PHARE PHARE-threads
    git -C PHARE-threads worktree add ../PHARE-master dc62fdd4e
    git -C PHARE-master submodule update --init --recursive    # if subprojects/ is empty

Use the way PHARE is usually built on GENOA in this account. What matters:

- `-DCMAKE_BUILD_TYPE=Release`, same flags for both builds. Build on (or for) a GENOA node.
- Only the 2D, first order module is needed: write `2,1,4` in a file and pass
  `-DPHARE_PERMUTATIONS=/path/to/that/file`. Targets: `cpp_2_1_4 cpp_etc dictator`.
- SAMRAI in Release too. A SAMRAI built as Debug or RelWithDebInfo makes every exchange call
  MPI barriers and spoils the timings. Part A below checks this.
- HDF5 with MPI support, as for any PHARE run (parts G and H write files).
- Python with numpy and h5py, and whatever `pyphare` needs (`requirements.txt`).

Then copy `env.sh.template` to `env.sh` and fill it in.

## Step 2: check the launch line before anything long

`run.sh` launches every run with

    srun --ntasks=R --cpus-per-task=T --hint=nomultithread --cpu-bind=cores --distribution=block:block

What the runs need, whatever the flags:

1. One hardware thread per core: a process must not be given both hardware threads of a core.
2. Each process is bound to its own T cores and to no other. PHARE then pins one thread on each
   of them. If a process may run on all the cpus of the node, PHARE does not pin anything and
   the timings are meaningless.
3. When R x T = 192 ("full"): process r has cores r*T to (r+1)*T - 1, so 8 x 24 is one process
   per group of 24 cores, 2 x 96 one process per socket.
4. When R x T < 192 ("spread", used by the scaling parts): the processes are shared out evenly
   between the two sockets.

Submit parts A, B and H first (about 15 minutes):

    sbatch --export=ALL,PARTS="A B H" job.sbatch      # after setting the account in job.sbatch

Then look at:

- `$RESULTS/affinity.txt`: one line per process. For the branch with T threads, `pinned=T` and
  `pinned_cpus` must be T different cores of the process, different from one process to the
  next. With 1 thread nothing is pinned by PHARE: `others_may_run_on` must be a single core
  (one or two cpu numbers, never a long range).
- `$RESULTS/machine.txt`: under "SAMRAI checks compiled in" nothing must be printed.
- `$RESULTS/bench.txt`: no `FAILED` line, and the three `CHECK` lines of each level as expected
  (see part H below).

If the binding is not as described, define `launch()` in `env.sh` (see the template) with the
flags that give it, run part B again, and report the launch line you ended with. The binding
wrapper CINES documents (`adastra_acc_binding.sh`, with `numactl --physcpubind`) is another way
to get it.

## Step 3: the measurements

    sbatch --export=ALL,PARTS="C" job.sbatch          # about 45 minutes, the main result
    sbatch --export=ALL,PARTS="D E F G" job.sbatch    # about 1 hour

The durations are estimates from a 48-core machine. Every run writes its line as soon as it
ends, so a job that hits its time limit loses one run only: submit the missing part again.

| part | what it measures | runs |
|------|------------------|------|
| A | the machine and the builds: `lscpu`, memory domains, modules, compiler, commits | none |
| B | where the threads run, for every layout used below | 14 short ones |
| H | same numbers as master? output files compared on 2 processes | 8 small ones |
| C | **main table**: master on 192 processes against the branch on 192x1, 96x2, 48x4, 32x6, 24x8, 16x12, 8x24, 4x48, 2x96, 1x192; uniform mesh (60 steps) and one refined level (15 steps); twice | 44 |
| D | master alone, 192 down to 1 process (20 steps): how MPI scales on this node | 9 |
| E | branch, 8 processes with 1 to 24 threads, 2 processes with 12 to 96 threads (20 steps) | 10 |
| F | patch size 25, 100, 200, 400 (C has 50): master, and branch on 16x12, 8x24, 2x96 | 16 |
| G | closer to production: tagged refinement regridded as the run goes, fields written every 10 steps: master, branch on 192x1, 32x6, 16x12, 8x24; twice | 10 |

Part H, expected lines (they were obtained on a workstation with 2 processes):

- master vs branch with 1 thread: `different 0` (the branch with one thread is master's code path);
- branch with 4 vs 12 threads: `different 0` (results do not depend on the number of threads);
- master vs branch with 12 threads: a worst relative difference of about 1e-15 (sums done in
  another order).

Two processes only, because master itself does not give the same last digits twice with three
or more.

## What to send back

The directory `$RESULTS` without `out/` (text only, well under 1 MB):

- `bench.txt`: one `BENCH ...` line per run, the `CHECK` and `FAILED` lines;
- `affinity.txt`, `machine.txt`, `logs/`, the `phare-threads-*.out` files of the jobs;
- `env.sh` and `job.sbatch` as used, and a few lines on what had to be changed and why.

Things worth a sentence in the report if you see them:

- a layout whose two repeats differ by more than 5 %;
- a run much slower than its neighbours (192x1 on the branch should be close to master);
- anything `machine.txt` shows about how a socket is split in memory domains;
- warnings of MPI about the thread level, or about binding.

## Not asked for now

Timings of each part of a time step (particles, fields, exchanges) need a separate build with
probes. They will be asked for in a second round, once these totals show where to look.
