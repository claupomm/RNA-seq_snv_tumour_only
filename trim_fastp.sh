#!/bin/bash

#SBATCH -c 4          # number of core to be used
#SBATCH -t 0-01:00    # estimated run-time in D-HH:MM
#SBATCH -p normal     # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=10G     # Memory pool for all cores

# Get sample name
sample=${PWD##*/}


# trim sequences
# -c, --correction                     enable base correction in overlapped regions (only for PE data), default is disabled
# -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
date
echo "Trim sequences for $sample..."
fastp -i *_1.f*q.gz -I *_2.f*q.gz -o ../../trim/$sample.R1.fq.gz -O ../../trim/$sample.R2.fq.gz -q 20 -c -p --thread 4
echo "Done."

mv fastp.html ../../qc/$sample.fastp.html
mv fastp.json ../../qc/$sample.fastp.json


