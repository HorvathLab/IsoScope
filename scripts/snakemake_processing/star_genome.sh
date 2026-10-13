#!/bin/bash

#SBATCH --job-name=star_genome
#SBATCH --time=01:30:00
#SBATCH -p cpu
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/star_genome2/%x_%j.out
#SBATCH --error=logs/star_genome2/%x_%j.err


module load star/2.7.11b

STAR --runThreadN "$1" \
        --runMode genomeGenerate \
        --genomeDir "$2" \
        --genomeFastaFiles "$3" \
        --sjdbGTFfile "$4" \
        --sjdbOverhang "$5"
