#!/bin/bash
#SBATCH -J ascot5      #Job name
#SBATCH -o %x.o%j        #stdout (%x=jobname, %j=jobid)
#SBATCH -e %x.e%j        #stderr (%x=jobname, %j=jobid)

#SBATCH --partition=standard     #Queue/Partition
#SBATCH --nodes=8            #Total number of nodes
#SBATCH --ntasks-per-node=1      #MPI tasks per node
#SBATCH --cpus-per-task=16	#CPUs per task for OpenMP
#SBATCH --time=2:59:00       #Wall clock limit
#SBATCH --mem=5GB                #Set mem./node requirement (default: 63000 MB, max: 190GB)

#SBATCH --mail-type=all       #Send mail, e.g. for begin/end/fail/none
#SBATCH --mail-user=linvelgal@gmail.com   #Mail address

echo Job name $SLURM_JOB_NAME
echo Job ID $SLURM_JOB_id

export OMP_NUM_THREADS=48

# Run the program:   
srun ./ascot5_mpi  --in $1
date
