#!/bin/bash

#SBATCH -c 14           # number of core to be used
#SBATCH -t 0-01:00      # estimated run-time in D-HH:MM
#SBATCH -p normal       # p=short <6h, p=mid <2d, p=long <4d
#SBATCH --mem=35G      	# Memory pool for all cores 


# Get sample name
sample=${PWD##*/}

# get paths, tools
star=/path/to/STAR
genome=/path/to/your/star/indices # genome with chr and scaffolds, primary assembly, transcriptome as recommended in Star manual for 2.7.10b

date
echo "Alignment in 2passmode for gatk rna-seq snv analysis..."
$star \
--genomeDir $genome \
--runThreadN 14 \
--readFilesIn ../../trim/$sample.R*.gz \
--readFilesCommand "gunzip -c" \
--sjdbOverhang 149 \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--limitBAMsortRAM 35000000000 \
--outSAMunmapped Within \
--outFileNamePrefix ../../star/$sample.
echo "Done."
date

mv ../../star/$sample.*bam ../../bam/$sample.bam
# samtools index ../../bam/$sample.bam

rm -r ../../star/$sample._STAR*
