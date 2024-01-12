#!/bin/bash

#SBATCH -c 8           # number of core to be used
#SBATCH -t 0-08:00     # estimated run-time in D-HH:MM
#SBATCH --mem=15G      # Memory pool for all cores




# Get sample name
sample=${PWD##*/}
gatk=~/Programme/gatk/gatk # -4.3.0.0
fasta=../data/Homo_sapiens_assembly38.fasta
indel=../data/Mills_and_1000G_gold_standard.indels.hg38.vcf
snp=../data/1000G_phase1.snps.high_confidence.hg38.vcf
chr_list=../data/chrName_head25.list


cd ../../bam


# Main Steps
# Data Cleanup
# Tools involved:
# MarkDuplicates
# We use MarkDuplicates (similarly to our DNA pre-processing best practices pipeline
# Preprocessing
# insert read group required for further processing
date
echo "Read Group adding..."
$gatk --java-options "-XX:ParallelGCThreads=8 \
	-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx10g" \
AddOrReplaceReadGroups \
    -I $sample.bam \
    -O $sample.rg.bam \
    -SORT_ORDER coordinate \
    -RGID bc.eurofins.$sample \
    -RGLB barstrand_specific \
    -RGPL illumina \
    -RGSM $sample \
    -RGPU flowcell.$sample \
    -CREATE_INDEX True \
    --verbosity ERROR
echo "Done."
echo ""

# MarkDuplicates run on a local Spark cluster with 8 executor cores
date
echo "Mark Duplicates..."
$gatk --java-options "-Xms10g" \
MarkDuplicatesSpark \
    -I $sample.rg.bam \
    -O $sample.mark_dup.bam \
    -M $sample.mark_dup_metrics.txt \
    -L $chr_list \
    -- \
    --spark-runner LOCAL \
    --spark-master 'local[8]' \
    --verbosity ERROR
echo "Done."
echo ""
rm $sample.rg.ba*

# SplitNCigarReads
date
echo "SplitNCigarReads..."
$gatk --java-options "-XX:ParallelGCThreads=8 \
	-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx10g" \
	SplitNCigarReads \
	-R $fasta \
	-L $chr_list \
	-I $sample.mark_dup.bam \
	-O $sample.cigar.bam \
	--verbosity ERROR
echo "Done."
echo ""
rm $sample.mark_dup.ba*

# Base Quality Recalibration
date
echo "Base (Quality Score) Recalibration..."
$gatk --java-options "-XX:ParallelGCThreads=8 \
	-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx10g" \
	BaseRecalibrator \
	-R $fasta \
	-I $sample.cigar.bam  \
	--use-original-qualities \
	--known-sites $snp \
	--known-sites $indel \
	-L $chr_list \
	-O $sample.recal.table \
	--verbosity ERROR
echo "Done."
echo ""

date
echo "Apply bqsr..."
$gatk --java-options "-XX:ParallelGCThreads=8 \
	-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx10g" \
	ApplyBQSR \
	--add-output-sam-program-record \
	-R $fasta \
	-I $sample.cigar.bam  \
	--use-original-qualities \
	-O $sample.recal.bam \
	--bqsr-recal-file $sample.recal.table \
	--verbosity ERROR
echo "Done."
echo ""
rm $sample.cigar.ba*
samtools index $sample.recal.bam

date



