#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J ComputeGroundMotionParameters
#Output and error (also --output, --error):
#SBATCH -o ./%j.%x.out
#SBATCH -e ./%j.%x.err
#SBATCH --chdir=***your working dir***
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=***your mail***
# Wall clock limit:
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=ALL
#SBATCH --account=***your account***
#constraints are optional
#--constraint="scratch&work"
#SBATCH --partition=test
#Number of nodes and MPI tasks per node: 
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --ear=off
#Run the program: 
export MP_SINGLE_THREAD=no
export OMP_NUM_THREADS=48
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS


myfile=***path to a -surface.xdmf file***
srun python -u SeisSol/postprocessing/science/GroundMotionParametersMaps/ComputeGroundMotionParametersFromSurfaceOutput_Hybrid.py $myfile --MP 48
