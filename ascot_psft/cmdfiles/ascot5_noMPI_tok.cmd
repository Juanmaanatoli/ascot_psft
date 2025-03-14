#!/bin/bash
#SBATCH -J ascot5      #Job name
#SBATCH -o ./%x.o%j        #stdout (%x=jobname, %j=jobid)
#SBATCH -e ./%x.e%j        #stderr (%x=jobname, %j=jobid)

#SBATCH --partition=p.tok.openmp     #Queue/Partition
#SBATCH --qos=p.tok.openmp.48h   #Quality of Service (see below): s.tok.short, s.tok.standard, s.tok.long, tok.debug
#SBATCH --nodes=1            #Total number of nodes
#SBATCH --ntasks-per-node=1      #MPI tasks per node
#SBATCH --cpus-per-task=32	#CPUs per task for OpenMP
#SBATCH --time=47:59:00       #Wall clock limit
#SBATCH --mem 10GB                #Set mem./node requirement (default: 63000 MB, max: 190GB)
##
#SBATCH --mail-type=all       #Send mail, e.g. for begin/end/fail/none
#SBATCH --mail-user=linvelgal@gmail.com   #Mail address

echo Job name $SLURM_JOB_NAME
echo Job ID $SLURM_JOB_id


# Run the program:  
date
srun ./ascot5_noMPI  --in $1 
date
