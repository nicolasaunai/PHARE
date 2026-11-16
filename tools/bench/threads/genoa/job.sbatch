#!/bin/bash
#SBATCH --job-name=phare-threads
#SBATCH --account=CHANGE_ME
#SBATCH --constraint=GENOA
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=02:00:00
#SBATCH --output=phare-threads-%j.out

# sbatch --export=ALL,PARTS="A B H C" job.sbatch      (from this directory)
cd "$SLURM_SUBMIT_DIR" && bash ./parts.sh ${PARTS:-A B H}
