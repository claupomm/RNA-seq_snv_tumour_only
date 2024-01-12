#!/bin/bash

#SBATCH -c 4             # number of core to be used
#SBATCH -t 0-02:00       # estimated run-time in D-HH:MM
#SBATCH --mem=15G        # Memory pool for all cores


# Get sample name
sample=${PWD##*/}
# variables
fasta=../data/Homo_sapiens_assembly38.fasta

cd ../../snp


# + convert vcf to maf for genvisr
# vep: ~/anaconda3/envs/py3.7/bin/ => ensembl-vep-105 => gnomad r2.1.1 >125k exomes
# vcf2maf.pl => default max_subpop_af > 0.0004 any gnomad population
# here: --max-subpop-af 0.01 => af > 0.01


# no filter (rna-edit, lcr sites)
vcf=$sample.hc.pass.vcf.gz
tmp=${vcf/.gz/}
gunzip -c $vcf > $tmp
perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $tmp --output-maf $tmp.maf --ref-fasta $fasta --tumor-id $sample --max-subpop-af 0.01 --ncbi-build GRCh38 --vep-path ~/anaconda3/envs/py3.7/bin/
pigz -p 4 $tmp.maf
pigz -p 4 ${tmp/vcf/vep.vcf} 
rm $tmp


# with lcr, without rna-edit sites
vcf=$sample.hc.pass2.vcf.gz
tmp=${vcf/.gz/}
gunzip -c $vcf > $tmp
perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $tmp --output-maf $tmp.maf --ref-fasta $fasta --tumor-id $sample --ncbi-build GRCh38 --max-subpop-af 0.01 --vep-path ~/anaconda3/envs/py3.7/bin/
pigz -p 4 $tmp.maf
pigz -p 4 ${tmp/vcf/vep.vcf} 
rm $tmp


# without lcr, without rna-edit sites
vcf=$sample.hc.pass2.lcr.vcf.gz
tmp=${vcf/.gz/}
gunzip -c $vcf > $tmp
perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $tmp --output-maf $tmp.maf --ref-fasta $fasta --tumor-id $sample --ncbi-build GRCh38 --max-subpop-af 0.01 --vep-path ~/anaconda3/envs/py3.7/bin/
pigz -p 4 $tmp.maf
pigz -p 4 ${tmp/vcf/vep.vcf} 
rm $tmp


# without lcr, without rna-edit sites, gatk 1000 genomes
vcf=$sample.hc.pass2.lcr.1kG.vcf.gz
tmp=${vcf/.gz/}
gunzip -c $vcf > $tmp
perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $tmp --output-maf $tmp.maf --ref-fasta $fasta --tumor-id $sample --ncbi-build GRCh38 --max-subpop-af 0.01 --vep-path ~/anaconda3/envs/py3.7/bin/
pigz -p 4 $tmp.maf
pigz -p 4 ${tmp/vcf/vep.vcf} 
rm $tmp


# without lcr, without rna-edit sites, nih dbsnp
vcf=$sample.hc.pass2.lcr.dbsnp.vcf.gz
tmp=${vcf/.gz/}
gunzip -c $vcf > $tmp
perl ~/Programme/mskcc-vcf2maf-754d68a/vcf2maf.pl --input-vcf $tmp --output-maf $tmp.maf --ref-fasta $fasta --tumor-id $sample --ncbi-build GRCh38 --max-subpop-af 0.01 --vep-path ~/anaconda3/envs/py3.7/bin/
pigz -p 4 $tmp.maf
pigz -p 4 ${tmp/vcf/vep.vcf} 
rm $tmp
