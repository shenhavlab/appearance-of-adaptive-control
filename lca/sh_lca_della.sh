#!/bin/bash

#SBATCH -t 12:00:00
#SBATCH --mem=4GB
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-type END
#SBATCH --mail-user hritz@princeton.edu
#SBATCH -J LCA-fit
#SBATCH --output slurm-logs/%A-%a_lca_fit.txt 
#SBATCH --array=1-5

module load matlab/R2024b


declare -A labels=(\
                    [1]="0.5"\
                    [2]="1.0"\
                    [3]="2.0"\
                    [4]="2.5"\
                    [5]="5"\
                    )


max_time=${labels[${SLURM_ARRAY_TASK_ID}]}
echo "max time: ${max_time}"

echo 'started at:'
date


matlab -nodisplay -nosplash -r "della_ibs(${max_time})" 


echo 'completed at:'
date