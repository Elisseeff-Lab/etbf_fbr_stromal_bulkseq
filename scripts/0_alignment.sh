#!/bin/bash

################################################################################
# Bulk RNA sequencing alignment pipeline
#
# This script performs alignment of raw FASTQ files to a reference genome.
#
# Data:
# - Raw FASTQ files are available on GEO under accession: GSEXXXXX
# - The FASTQ files are gzipped and named in the format: <sample_name>_R1.fastq.gz and <sample_name>_R2.fastq.gz
# - The reference genome uses the GRCm39 primary assembly and the GENCODE M27
# annotation. They can be downloaded here: https://www.gencodegenes.org/mouse/release_M27.html
#
# Modify these paths before running:
# WORKING_DIR: Working directory; it is assumed FASTQ files are in a subdirectory called "fastqs"
# REFERENCE_DIR: Directory containing the reference genome indexes
# 
# Output:
# - The script generates STAR output files for each sample, including:
#   - SAM files
#   - Gene counts (used for downstream analysis)
#   - Log files
################################################################################

# Load modules
module load star # STAR v2.7.8a used; change module or activation as needed

WORKING_DIR="/path/to/working/directory" # Modify this path
REFERENCE_DIR="/path/to/reference/directory" # Modify this path

# Grab sample names by splitting file name at underscore
cd ${WORKING_DIR}/fastqs
SAMPLES=()
for F in *
do
    VAR="$(cut -d'_' -f1 <<<$F)"
    SAMPLES+=($VAR)
done

# Remove duplicates (since we have paired end samples)
SAMPLES=($(echo "${SAMPLES[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

cd ${WORKING_DIR}

# Align all samples using STAR
for SAMP in ${SAMPLES[@]}
do
echo "Working on $SAMP"
STAR --runThreadN 8 \
--genomeDir ${REFERENCE_DIR} \
--readFilesCommand zcat \
--readFilesIn fastqs/${SAMP}_R1.fastq.gz fastqs/${SAMP}_R2.fastq.gz \
--quantMode GeneCounts \
--outFileNamePrefix ${SAMP}
done

