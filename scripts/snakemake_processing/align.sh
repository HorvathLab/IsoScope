#!/bin/bash

#SBATCH --job-name=star
#SBATCH --time=12:00:00
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --output=logs/star_short2/%x_%j.out
#SBATCH --error=logs/star_short2/%x_%j.err

module load star

STAR --runThreadN "$1" \
        --genomeDir "$2" \
        --readFilesIn "$3" "$4" \
        --twopassMode Basic \
        --sjdbGTFfile /SMHS/groups/mgpc/horvathlab/ml/rerun/new_refs/gencode.v49.basic.annotation.gtf \
        --outSAMtype BAM Unsorted \
        --quantMode GeneCounts TranscriptomeSAM \
        --outSAMattrRGline ID:1 LB:lib PL:illumina PU:unit SM:smp \
        --outFileNamePrefix "$5" \
        --outSAMmapqUnique 60
