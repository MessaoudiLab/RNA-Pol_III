#!/bin/bash
  
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=15:00:00    # 1 day and 15 minutes
#SBATCH --mail-user=elizabethvance03@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="FeatureCounts"

# $1 = path to aligned and sorted bam files
# $2 = path to annotation gtf file
# $3 = name of output file

fileList=$(ls -d "$1"*.bam)
singularity exec ./singularity/rnaSeq.sif bash featureCounts.sh $fileList $2 $3
