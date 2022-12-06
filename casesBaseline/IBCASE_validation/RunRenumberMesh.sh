#!/bin/bash
#SBATCH -D .  # Working directory
#SBATCH --job-name=sprayC_BL_corr                # Job name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=raffaele.bellini0@city.ac.uk         # Where to send mail	
#SBATCH --exclusive                          # Exclusive use of nodes
#SBATCH --nodes=3                            # Run on 2 nodes (each node has 48 cores)
#SBATCH --ntasks-per-node=48                 # Use all the cores on each node
#SBATCH --mem=0                              # Expected memory usage (0 means use all available memory)
#SBATCH --time=2-23:59:00                      # Time limit hrs:min:sec
#SBATCH --output=cluster100.out        # Standard output and error log [%j is replaced with the jobid]
#SBATCH --error=cluster100.error

#enable modules
source /opt/flight/etc/setup.sh
flight env activate gridware

srun hostname  | sort > hosts.$SLURM_JOB_ID

#Command line to run task
date > dateLog
mpirun -np 144 renumberMesh -overwrite  -parallel  1>log_renumber>&1


date > dateLog

