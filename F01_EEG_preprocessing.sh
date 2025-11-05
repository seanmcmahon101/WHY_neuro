#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --time 360:0
#SBATCH --mem 50G
#SBATCH --qos bbdefault
#SBATCH --mail-type ALL
#SBATCH -o logs/EEG_preprocessing-%j.out
#SBATCH --array=1-36

set -eu

module purge
module load bluebear
module load MATLAB/2020a

matlab -nodisplay -r "run F01_preprocessing(${SLURM_ARRAY_TASK_ID}), quit"

