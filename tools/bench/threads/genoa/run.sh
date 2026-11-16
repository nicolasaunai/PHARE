#!/bin/bash
# One run of harris_2d.py inside a SLURM allocation of one node:
#
#   run.sh <master|branch> <ranks> <threads per rank> <label> [full|spread]
#
# The problem is set by the BENCH_* variables of harris_2d.py. The BENCH line of the run is
# appended to $RESULTS/bench.txt with the label and the build, the whole output is kept in
# $RESULTS/logs/. A run that fails or times out adds a FAILED line and the sequence goes on.
#
# full   : ranks x threads fills the node, rank r gets cores [r*threads, (r+1)*threads)
# spread : fewer cores than the node has, ranks are shared out between the two sockets

here=$(dirname "$(readlink -f "$0")")
source "$here/env.sh"

build=$1 ranks=$2 threads=$3 label=$4 placement=${5:-full}

case $build in
master) src=$PHARE_MASTER_SRC bld=$PHARE_MASTER_BUILD ;;
branch) src=$PHARE_BRANCH_SRC bld=$PHARE_BRANCH_BUILD ;;
*) echo "run.sh: build is master or branch" && exit 2 ;;
esac

# one hardware thread per core; a rank is bound to its own cores, that it may not leave
if ! declare -F launch > /dev/null; then
    launch() {
        local distribution=block:block
        [ "$3" = spread ] && distribution=block:cyclic
        srun --ntasks="$1" --cpus-per-task="$2" --hint=nomultithread --cpu-bind=cores \
            --distribution=$distribution "${@:4}"
    }
fi

export PYTHONPATH=$bld:$src/pyphare:$src
export PHARE_THREADS=$threads # read by the branch only
export OMP_NUM_THREADS=1
export BENCH_TAG="label=$label build=$build placement=$placement"

mkdir -p "$RESULTS/logs"
n=$(find "$RESULTS/logs" -name '*.log' | wc -l)
log=$RESULTS/logs/$(printf %03d "$n")-$label-$build-${ranks}x$threads.log

echo "RUN $BENCH_TAG ranks=$ranks threads=$threads $(env | grep ^BENCH_ | grep -v BENCH_TAG | sort | tr '\n' ' ')" > "$log"
# the same script for both builds: master has no tools/bench/threads
export -f launch
timeout "${RUN_TIMEOUT:-1800}" bash -c 'launch "$@"' launch "$ranks" "$threads" "$placement" \
    python3 -u "$PHARE_BRANCH_SRC/tools/bench/threads/harris_2d.py" >> "$log" 2>&1
status=$?

if [ $status -eq 0 ] && grep -q "^BENCH " "$log"; then
    grep "^BENCH " "$log" | tee -a "$RESULTS/bench.txt"
    grep "^BENCH_AFFINITY" "$log" | sort -t= -k2 -n | sed "s/^BENCH_AFFINITY/BENCH_AFFINITY $BENCH_TAG ranks=$ranks/" >> "$RESULTS/affinity.txt"
else
    echo "FAILED $BENCH_TAG ranks=$ranks threads=$threads exit=$status log=$(basename "$log")" | tee -a "$RESULTS/bench.txt"
fi
