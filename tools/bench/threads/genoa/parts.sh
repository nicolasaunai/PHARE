#!/bin/bash
# The sequence of runs, by parts, inside a SLURM allocation of one whole node:
#
#   parts.sh A B H C        (see README.md for what each part measures)
#
# Every run appends one line to $RESULTS/bench.txt as soon as it ends: a job that runs out of
# time loses the run it was in, nothing else. A part can be run again, its lines are appended.

here=$(dirname "$(readlink -f "$0")")
source "$here/env.sh"
run=$here/run.sh
mkdir -p "$RESULTS"

CORES=${CORES:-192}
# threads per rank of the layouts that fill the node: 192x1, 96x2, ... 2x96, 1x192
THREADS=${THREADS:-"1 2 4 6 8 12 24 48 96 192"}
REPEATS=${REPEATS:-2}

# the problem: 2D Harris sheet, 128 million particles, 512 patches on the coarse level
BIG="BENCH_CELLS=${BIG_CELLS:-1600x800} BENCH_PPC=${BIG_PPC:-100} BENCH_PATCH=${BIG_PATCH:-50}"
UNIFORM="$BIG BENCH_LEVELS=1 BENCH_STEPS=${UNIFORM_STEPS:-60}"
REFINED="$BIG BENCH_LEVELS=2 BENCH_STEPS=${REFINED_STEPS:-15}"
SCALING="$BIG BENCH_LEVELS=1 BENCH_STEPS=${SCALING_STEPS:-20}"
# small ones: to look at the binding, and to compare output files
TINY="BENCH_CELLS=${TINY_CELLS:-800x400} BENCH_PPC=5 BENCH_PATCH=25 BENCH_LEVELS=1 BENCH_STEPS=2"
CHECK="BENCH_CELLS=${CHECK_CELLS:-400x200} BENCH_PPC=20 BENCH_PATCH=50 BENCH_STEPS=10"

part_A() { # the machine and the two builds
    {
        echo "== date"; date
        echo "== lscpu"; lscpu
        echo "== memory domains (cpus of each)"; grep . /sys/devices/system/node/node*/cpulist
        echo "== numactl"; numactl -H 2>&1
        echo "== memory"; free -g
        echo "== modules"; module list 2>&1
        echo "== srun"; srun --version
        echo "== python"; python3 --version; python3 -c "import numpy, h5py; print(numpy.__version__, h5py.__version__)"
        for build in BRANCH MASTER; do
            src=PHARE_${build}_SRC bld=PHARE_${build}_BUILD
            echo "== $build source"; git -C "${!src}" log --oneline -1; git -C "${!src}" status --short | head
            echo "== $build build"
            grep -E "^(CMAKE_BUILD_TYPE|CMAKE_CXX_COMPILER|CMAKE_CXX_FLAGS|CMAKE_CXX_FLAGS_RELEASE|withIPO|devMode|SAMRAI_ROOT|MPI_CXX_COMPILER)[:=]" "${!bld}/CMakeCache.txt"
            echo "== $build SAMRAI checks compiled in (must print nothing)"
            find "${!bld}" $(grep -E "^SAMRAI_ROOT" "${!bld}/CMakeCache.txt" | cut -d= -f2) -name SAMRAI_config.h 2> /dev/null |
                xargs -r grep -H "define DEBUG_CHECK_ASSERTIONS"
            echo "== $build libraries of the module"; ldd "${!bld}"/pybindlibs/cpp_2_1_4*.so | grep -E -i "mpi|hdf5|SAMRAI_tbox|stdc"
        done
    } > "$RESULTS/machine.txt" 2>&1
    echo "part A: $RESULTS/machine.txt"
}

part_B() { # where the threads run, for every layout used below: a few seconds each
    for t in $THREADS; do
        r=$((CORES / t))
        [ "$t" = 1 ] && env $TINY BENCH_AFFINITY=1 $run master "$r" 1 binding
        env $TINY BENCH_AFFINITY=1 $run branch "$r" "$t" binding
    done
    env $TINY BENCH_AFFINITY=1 $run master 8 1 binding spread
    env $TINY BENCH_AFFINITY=1 $run branch 8 6 binding spread
    env $TINY BENCH_AFFINITY=1 $run branch 2 12 binding spread
}

part_C() { # the main table: every layout that fills the node, against master on all cores
    for repeat in $(seq "$REPEATS"); do
        for problem in uniform refined; do
            settings=$UNIFORM && [ $problem = refined ] && settings=$REFINED
            env $settings $run master "$CORES" 1 $problem
            for t in $THREADS; do
                env $settings $run branch $((CORES / t)) "$t" $problem
            done
        done
    done
}

part_D() { # how master scales with ranks alone, same problem
    for r in $CORES $((CORES / 2)) $((CORES / 4)) $((CORES / 6)) $((CORES / 12)) $((CORES / 24)) 4 2 1; do
        [ "$r" -ge 1 ] || continue
        placement=spread && [ "$r" = "$CORES" ] && placement=full
        env $SCALING RUN_TIMEOUT=3000 $run master "$r" 1 mpi-scaling $placement
    done
}

part_E() { # how the branch scales with threads, for 8 ranks and for 2 ranks
    for t in 1 2 4 8 12 24; do
        placement=spread && [ $((8 * t)) = "$CORES" ] && placement=full
        env $SCALING $run branch 8 "$t" thread-scaling $placement
    done
    for t in 12 24 48 96; do
        placement=spread && [ $((2 * t)) = "$CORES" ] && placement=full
        env $SCALING $run branch 2 "$t" thread-scaling $placement
    done
}

part_F() { # patch size: few large patches are the weak case of threads over patches
    for patch in 25 100 200 400; do
        env $UNIFORM BENCH_PATCH=$patch $run master "$CORES" 1 patch-size
        for t in 12 24 96; do
            env $UNIFORM BENCH_PATCH=$patch $run branch $((CORES / t)) "$t" patch-size
        done
    done
}

part_G() { # closer to a production run: cells tagged, hierarchy regridded, fields written
    settings="$BIG BENCH_LEVELS=2 BENCH_REFINE=tagging BENCH_STEPS=${REGRID_STEPS:-30} BENCH_DIAG_EVERY=10"
    for repeat in $(seq "$REPEATS"); do
        env $settings BENCH_DIAG_DIR=$RESULTS/out/regrid $run master "$CORES" 1 regrid-output
        for t in 1 6 12 24; do
            env $settings BENCH_DIAG_DIR=$RESULTS/out/regrid $run branch $((CORES / t)) "$t" regrid-output
        done
    done
    rm -rf "$RESULTS/out/regrid"
}

part_H() { # same numbers as master? 2 ranks: master itself is not reproducible with 3 or more
    out=$RESULTS/out
    for levels in 1 2; do
        env $CHECK BENCH_LEVELS=$levels BENCH_DIAG_DIR=$out/master $run master 2 1 check spread
        env $CHECK BENCH_LEVELS=$levels BENCH_DIAG_DIR=$out/branch1 $run branch 2 1 check spread
        env $CHECK BENCH_LEVELS=$levels BENCH_DIAG_DIR=$out/branch4 $run branch 2 4 check spread
        env $CHECK BENCH_LEVELS=$levels BENCH_DIAG_DIR=$out/branch12 $run branch 2 12 check spread
        {
            echo "CHECK levels=$levels master vs branch, 1 thread (expected: different 0): $(python3 "$here/compare.py" $out/master $out/branch1 | tail -1)"
            echo "CHECK levels=$levels branch, 4 vs 12 threads (expected: different 0): $(python3 "$here/compare.py" $out/branch4 $out/branch12 | tail -1)"
            echo "CHECK levels=$levels master vs branch, 12 threads (expected: about 1e-15): $(python3 "$here/compare.py" $out/master $out/branch12 | tail -1)"
        } | tee -a "$RESULTS/bench.txt"
        rm -rf "$out"
    done
}

for part in "$@"; do
    echo "PART $part $(date +%H:%M:%S)" | tee -a "$RESULTS/bench.txt"
    part_"$part"
done
echo "DONE $* $(date +%H:%M:%S)" | tee -a "$RESULTS/bench.txt"
