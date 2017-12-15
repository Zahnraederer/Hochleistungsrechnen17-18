#!/bin/bash
	
# Time limit is one minute. See "man sbatch" for other time formats.
#SBATCH --time=1
# Run a total of ten tasks on two nodes (that is, five tasks per node).
#SBATCH --nodes=3
#SBATCH --ntasks=3
# Use "west" partition.
#SBATCH --partition=west
# Output goes to "job.out", error messages to "job.err".
#SBATCH --output=job.out
#SBATCH --error=job.err

. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" = 'abu' ]
then
	export MPICH_NEMESIS_NETMOD=ib
fi

#srun hostname
mpiexec ./partdiff-par 1 2 0 2 2 1
