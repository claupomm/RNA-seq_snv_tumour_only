#!/bin/bash

#SBATCH -c 4            # number of core to be used
#SBATCH -t 0-01:00      # estimated run-time in D-HH:MM
#SBATCH --mem=5G        # Memory pool for all cores


# Get sample name
sample=${PWD##*/}
# variables
filter=../data/1000G_phase3_v4_20130502.sites.hg38.af01.vcf
filter2=../data/GCF_000001405.40.common.chr.tsv

cd ../../snp



# filter out common snps from called snvs, without lcr, without rna-edit sites
i=$sample.hc.pass2.lcr.vcf.gz

# without indels for mutsigs, gatk 1000 genomes
j=$sample.hc.pass2.lcr.1kG.snp.vcf.gz
vcftools --gzvcf $i --exclude-positions $filter --remove-indels --recode --recode-INFO-all --stdout | pigz -p4 > $j
# with indels for waterfall, gatk 1000 genomes
j=$sample.hc.pass2.lcr.1kG.vcf.gz
vcftools --gzvcf $i --exclude-positions $filter --recode --recode-INFO-all --stdout | pigz -p4 > $j

# without indels for waterfall, nih dbsnp
j=$sample.hc.pass2.lcr.dbsnp.snp.vcf.gz
vcftools --gzvcf $i --exclude-positions $filter2 --remove-indels --recode --recode-INFO-all --stdout | pigz -p4 > $j
# with indels for waterfall, nih dbsnp
j=$sample.hc.pass2.lcr.dbsnp.vcf.gz
vcftools --gzvcf $i --exclude-positions $filter2 --recode --recode-INFO-all --stdout | pigz -p4 > $j


