#!/bin/sh
#SBATCH --job-name=parallel_cdrds
#SBATCH --output=run_par.out
#SBATCH --error=run_par.err
#SBATCH --account=pi-ejadams
#SBATCH --partition=broadwl
#SBATCH --time=3:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=28
#SBATCH --exclusive

cd /project2/ejadams/boughter/cdrds
# A LOT OF THIS STOLEN FROM RCC PAGE for "RUN-PARALLEL"
# Load the default version of GNU parallel.
module load parallel
#module load Anaconda3
module load python/anaconda-2020.02

source activate my_root

# MAKE SURE THIS ACTUALLY MATCHES NTASKS*NODES
ntasks=112

# When running a large number of tasks simultaneously, it may be
# necessary to increase the user process limit.
ulimit -u 10000

# This specifies the options used to run srun. The "-N1 -n1" options are
# used to allocates a single core to each task.
srun="srun --exclusive -N1 -n1"

# This specifies the options used to run GNU parallel:
#
#   --delay of 0.2 prevents overloading the controlling node.
#
#   -j is the number of tasks run simultaneously.
#
#   The combination of --joblog and --resume create a task log that
#   can be used to monitor progress.
#
parallel="parallel --delay 0.2 -j $SLURM_NTASKS --joblog runtask.log --resume"

# Aaand run it.
#$parallel "$srun python midway_null_bootstrap.py $ntasks {1} > test.{1}" ::: {1..112}

$parallel "$srun python midway_bootstrap_null.py $ntasks {1} > test.{1}" ::: {1..112}

# These jobs should all run in parallel, but if for some reason we're short on CPUs
# then they'll run in serial. Shouldn't happen though
