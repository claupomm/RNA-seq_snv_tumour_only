#!/bin/bash

#SBATCH -c 2             # number of core to be used
#SBATCH -t 0-12:00       # estimated run-time in D-HH:MM
#SBATCH --mem=15G        # Memory pool for all cores



# Get sample name
sample=${PWD##*/}
# variables
gatk=~/Programme/gatk/gatk # -4.3.0.0
fasta=../data/Homo_sapiens_assembly38.fasta
snp=../data/1000G_phase3_v4_20130502.sites.hg38.vcf
chr_list=../data/chrName_head25.list


cd ../../bam


# Variant Calling
# Run Haplotypecaller
date
echo "Mutation calling for RNA-Seq data via HaplotypeCaller..."
$gatk --java-options "-XX:ParallelGCThreads=8 -XX:GCHeapFreeLimit=10 -Xms10g" \
HaplotypeCaller -R $fasta \
    -L $chr_list \
    -I $sample.recal.bam \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp $snp \
    --dont-use-soft-clipped-bases \
    -O ../snp/$sample.hc.vcf \
	--verbosity ERROR
echo "Done."
echo ""


# minimal depth >=5
date
echo "Filter variants..."
$gatk --java-options "-XX:ParallelGCThreads=8 -XX:GCHeapFreeLimit=10 -Xmx10g" \
VariantFiltration \
    -R $fasta \
	-V ../snp/$sample.hc.vcf \
	-L $chr_list \
    --window 35 \
    --cluster 3 \
    --filter-name "DP" \
	--filter "DP < 5" \
    -O ../snp/$sample.hc.filt.vcf \
	--verbosity ERROR
echo "Done."
echo ""

# output with pass variant only
date
echo "Pass variants only to vcf..."
bcftools filter -i 'FILTER="PASS"' -o ../snp/$sample.hc.pass.vcf ../snp/$sample.hc.filt.vcf
echo "Done."
# echo ""


rm ../snp/$sample.hc.filt*.vcf ../snp/$sample.hc.filt*.vcf.idx
pigz -p 8 ../snp/$sample.hc.pass.vcf
pigz -p 8 ../snp/$sample.hc.vcf


date

