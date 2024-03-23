#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J compact
#Output and error (also --output, --error):
#SBATCH -o ./%j.%x.out
#SBATCH -e ./%j.%x.out

#Initial working directory:
#SBATCH --chdir=./

#Notification and type
#SBATCH --mail-type=END
# Wall clock limit:
#SBATCH --time=00:30:00
#SBATCH --no-requeue

#Setup of execution environment
#SBATCH --export=ALL
#SBATCH --account=pn49ha
#SBATCH --partition=test

#Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#EAR may impact code performance
#SBATCH --ear=off
#SBATCH --mem=80G

set -euo pipefail
mkdir -p extracted_output

counter=0
total_params=$(ls output/dyn_*-surface.xdmf | wc -l)
total_params=$((total_params * 3))

for current_file in output/dyn_*-surface.xdmf; do
    counter=$((counter+1))
    echo "Processing file $counter of $total_params: $current_file"
    srun -n 1 -c 1 --exclusive --mem-per-cpu 8G seissol_output_extractor $current_file --time "i:" --variable u1 u2 u3 --add2prefix _disp &
    # Improved check: avoids unnecessary wait on the first iteration
    #if (( $counter % $SLURM_NTASKS == 0 )); then
    if (( $counter % 10 == 0 )); then
      echo "waiting, $counter"
      wait
    fi
done


for current_file in output/*-fault.xdmf; do
    counter=$((counter+1))
    echo "Processing file $counter of $total_params: $current_file"
    srun -n 1 -c 1 --exclusive --mem-per-cpu 8G seissol_output_extractor $current_file &
    if (( $counter % 10 == 0 )); then
      echo "waiting, $counter"
      wait
    fi
done

for current_file in output/dyn_*-energy.csv; do
    counter=$((counter+1))
    echo "Processing file $counter of $total_params: $current_file"
    srun -n 1 -c 1  --exclusive --mem-per-cpu 8G cp $current_file extracted_output &
    if (( $counter % 10 == 0 )); then
      echo "waiting, $counter"
      wait
    fi
done
wait
mv *_extracted* extracted_output 
mv *_disp* extracted_output
