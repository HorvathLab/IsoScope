#!/bin/bash

#SBATCH --job-name=salmon
#SBATCH --time=6:00:00
#SBATCH -p cpu
#SBATCH --mem=128G
#SBATCH --output=logs/salmon2/%x_%j.out
#SBATCH --error=logs/salmon2/%x_%j.err

module load salmon/1.10.0

salmon quant \
-t "$1" \
-l A \
-a "$2" \
-p "$3" \
-o "$4" \

echo "job %x_%j"
