#!/bin/bash

#SBATCH --job-name=download_sra
#SBATCH --time=2:50:00
# SBATCH -p nano
#SBATCH -p cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --output=logs/download2/%x_%j.out
#SBATCH --error=logs/download2/%x_%j.err

BASEDIR="/SMHS/groups/mgpc/horvathlab/ml"
SRA_DIR="$BASEDIR/sra"
FASTQ_DIR="$BASEDIR/fastq"
TMP_DIR="$BASEDIR/tmp"

module load sratoolkit

prefetch --max-size 200GB "$1" --output-directory "$SRA_DIR"

sleep 60

fasterq-dump "$SRA_DIR/$1" \
        --outdir "$FASTQ_DIR" \
        --temp "$TMP_DIR" \
        --threads $2 \
        --split-files
