#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=15:00:00    # 1 day and 15 minutes
#SBATCH --mail-user=elizabethvance03@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="DESeq2"
singularity exec ./singularity/rnaSeq.sif Rscript DESEQ2_Messaoudi.R $1 $2 $3
