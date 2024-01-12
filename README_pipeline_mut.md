# Variant calling pipeline on RNA-seq data for tumour-only breast cancer cell lines

Typically, variant calling on tumour material is performed on whole-exome or whole-genome DNA sequencing. This workflow, however, has been used for analysing RNA-seq data of 29 tumour breast cancer cell lines without matched-normal sample pairs for identifying single nuleotide variants (SNVs) and insertions / deletions (InDels). Since capturing germline variants produces a massive amount of mutations per sample, the filtering process is of special importance in order to retrieve an essence of potentially intriguing variants. 

## Outline

1. Trimming via fastp
2. Aligning via STAR in two-pass mode
3. Group addition, read duplicate removal, SplitNCigarReads, base recalibration, germline calling via HaplotypeCaller for RNA-seq, variant filtering via the rich GATK tool bundle
4. Filtering of 
   1. RNA edit sites
   2. low complexity regions (LCRs)
   3. sites of low depth (<5)
   4. common variants described in 
      1. dbSNP
      2. 1000 Genomes
      3. gnomAD
   5. non-coding regions
   6. variants occuring in >20% of the breast cancer cell lines

## Intial steps

Following tools / programs need to be installed and running:

- fastp for trimming the reads
- STAR aligner for mapping
- GATK bundle of the Broad institute: removal of read duplicates, variant calling, etc
- samtools
- snpeff/snpsift
- vcftools
- VEP, ensembl-vep-105 includes gnomad r2.1.1
- vcf2maf: https://github.com/mskcc/vcf2maf/blob/main/README.md
- R-packages GenVisR, MutationalPatterns

Following data need to be downloaded:

- for GATK:
  
  - sources: 
    
    - https://gatk.broadinstitute.org/hc/en-us/articles/360035890811
    - or look here: https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&supportedpurview=project&prefix=&forceOnObjectsSortingFiltering=false
  
  - files 
    
    - Homo_sapiens_assembly38.fasta
    - 1000G_phase1.snps.high_confidence.hg38.vcf
    - Mills_and_1000G_gold_standard.indels.hg38.vcf

- RNA edit sites, http://srv00.recas.ba.infn.it/atlas/download.html: TABLE1_hg38_chr_pos.txt

- LCRs, https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true: LCR-hs38_with_chr.bed.gz

- COSMIC mutation data https://cancer.sanger.ac.uk/cell_lines/download:
  
  - CosmicCLP_MutantExport_v97_2023_01.tsv.gz
  - Cosmic_MutantCensus_v98_GRCh38.tsv.gz

- Common SNP sites
  
  - dbSNP, GCF_000001405.40.common.chr.tsv
  - 1000 Genomes, 1000G_phase3_v4_20130502.sites.hg38.af01.vcf
  - gnomAD: contained in VEP data sources



### Create folders

```
DIR=~/Documents/ngs/Project_bc_mut
mkdir $DIR
mkdir $DIR/raw
mkdir $DIR/data # store all downloaded data here
cd $DIR/raw
```



### Link sequencing files

Our RNA-seq fastq files are listed in file2sample.csv and can be downloaded from BioStudies: https://www.ebi.ac.uk/biostudies/studies/S-BSST1200

```
DIR_fastq=/path/to/fastq_files
for i in $(cut -d, -f1 ../file2sample.csv| tail -n+2 - ); do
if [[ $(ls $DIR_fastq) =~ ${i//-/_} ]]; then # check if sample is in subdirectory
mkdir $DIR/raw/$i
cd $DIR/raw/$i
ln -s $dir/*gz .
fi
done
```

## Preprocessing, mapping, and variant calling

Preprocessing steps best to be done on high computing cluster. Here, we do this via the job scheduler SLURM.

Start preprocessing:

```
cd $DIR
SAMPLES=$(find $DIR/raw/* -maxdepth 1 -type d)
mkdir star bam qc trim snp plots tables

# iterate each sample
for SAMPLE in $SAMPLES; do
cd $SAMPLE
sample=${PWD##*/}
# quality control and trimming via fastp
# for ~30M reads per sample RNA-seq < 10min
RES=$(sbatch -J $sample.1 -o $sample.1.out -e $sample.1.err ../../trim_fastp.sh)
# alignment via star, 2-pass mode <40min
RES2=$(sbatch --dependency=afterok:${RES##* } -J $sample.2 -o $sample.2.out -e $sample.2.err ../../align.sh)
# mark duplicates via picard + coverage, metrix, insert size, <4h
RES3=$(sbatch --dependency=afterok:${RES2##* } -J $sample.3 -o $sample.3.out -e $sample.3.err ../../gatk_preprocess.sh)
# HC for RNA-seq, several hours
sbatch --dependency=afterok:${RES3##* } -J $sample.4 -o $sample.4.out -e $sample.4.err ../../gatk_hc.sh
done
```

## Filtering steps

### Filter RNA-edit sites

Some basic knowledge about RNA-edit sites is explained here: https://academic.oup.com/nar/article/49/D1/D1012/5940507?login=true

```
rna_edit=~/Dokumente/ngs/Genomes/rediportal/TABLE1_hg38_chr_pos.txt
for SAMPLE in $SAMPLES; do
cd $SAMPLE
sample=${PWD##*/}
sbatch -J $sample.5 -o $sample.5.out -e $sample.5.err <<EOF
#!/bin/sh
cd ../../snp
vcftools --gzvcf $sample.pass.vcf.gz --recode --exclude-positions $rna_edit --stdout | gzip -c > $sample.pass2.vcf.gz
vcftools --gzvcf $sample.hc.pass.vcf.gz --recode --exclude-positions $rna_edit --stdout | gzip -c > $sample.hc.pass2.vcf.gz
EOF
done
```

### Filter low complexity regions (LCRs)

* A useful tutorial: https://github.com/hbctraining/variant_analysis/blob/main/lessons/10_variant_filtering.md

* Theory to LCRs: https://academic.oup.com/bioinformatics/article/30/20/2843/2422145?login=true
1. Get the bed file:

```
SNPEFF=/path/to/snpEff
cd $SNPEFF
curl -o LCR-hs38_with_chr.bed.gz -L https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true
gunzip -c LCR-hs38_with_chr.bed.gz > LCR-hs38_with_chr.bed
sed 's/^chr//g' LCR-hs38_with_chr.bed > LCR-hs38.bed
mv LCR-hs38.bed $DIR/data/.
```

2. Filter LCRs

```
cd $DIR/snp

for sample in $(cut -d, -f1 ../file2sample.csv| tail -n+2 - ); do
java -jar $SNPEFF/SnpSift.jar intervals \
-noLog -x -i $sample.hc.pass2.vcf.gz \
$SNPEFF/LCR-hs38.bed | gzip > $sample.hc.pass2.lcr.vcf.gz
done
```

### Filter common snps from called snvs

Prepare the dbSNP filter together with the downloaded dbsnp file at https://ftp.ncbi.nih.gov/snp/latest_release/VCF/ stored in $DIR/data

```
mkdir $DIR/data

# prepare file for filtering called variants
zgrep "^#" GCF_000001405.40.gz > GCF_000001405.40.common.vcf

# label COMMON for variants with AF>0.01 (1%)
zgrep ";COMMON" GCF_000001405.40.gz >> GCF_000001405.40.common.vcf

# get chr and position for common variants (>0.01)
refseq=$(cut -f1 gcf_refseq2chr.tsv)
IFS=' ' read -ra refseq <<< $(echo $refseq)
chr=$(cut -f2 gcf_refseq2chr.tsv)
IFS=$' ' read -ra chr <<< $(echo $chr)
touch GCF_000001405.40.common.chr.tsv

for i in ${!refseq[@]}; do
date
echo ${chr[$i]}
zgrep ${refseq[$i]} GCF_000001405.40.common.vcf.gz | sed "s/${refseq[$i]}/${chr[$i]}/g" | cut -f1,2 >> GCF_000001405.40.common.chr.tsv
done
```

Filter dbSNP and 1000 Genomes mutations with these data files:

- 1000G_phase3_v4_20130502.sites.hg38.af01.vcf
- GCF_000001405.40.common.chr.tsv

```
for sample in $(cut -d"," -f1 $DIR/file2sample.csv | tail -n29); do
cd $DIR/raw/$sample
sbatch -J $sample.6 -o $sample.6.out -e $sample.6.err ../../filter_common_variants.sh
done
```

### Convert vcf to maf

The maf file format is required for annotation, statistics, and visualisation. vcf2maf.pl will convert vcf to maf and filters gnomad with max af 0.01 for any population and includes VEP (ensembl-vep-105 => gnomad r2.1.1):

```
for sample in $(cut -d"," -f1 $DIR/file2sample.csv | tail -n29); do
cd $DIR/raw/$sample
sbatch -J $sample.7 -o $sample.7.out -e $sample.7.err ../../vcf2maf.sh
    done
```

### Annotate and summarise variants for all 29 breast cancer cell lines

```
R --file="mut_filter.R"
```

## Variant validation

### Statistics

Get numbers for filtered mutants for each step:

```
    cd $DIR/snp
    zgrep -c "^chr" *.hc.vcf.gz > stats_hc_all.txt # q20
    zgrep -c "^chr" *.hc.pass.vcf.gz > stats_hc_pass.txt # q20 + depth
    zgrep -c "^chr" *.hc.pass2.vcf.gz > stats_hc_pass2.txt # q20 + depth + rna_edit
    zgrep -c "^chr" *.hc.pass2.lcr.vcf.gz > stats_hc_lcr.txt # q20 + depth + rna_edit + lcr
    zgrep -c "^chr" *.hc.pass2.lcr.dbsnp.vcf.gz > stats_hc_dbsnp.txt # q20 + depth + rna_edit + lcr + dbsnp
```

Visualise variant number after each step of filtering:

```
    R --file="stats.R"
```

### Mutational signatures

Filter mutations by low complexity regions LCR, gnomad af > 0.01, extract single base substitutions (SBS). The R-package MutationalPatterns is used for calculating mutational signatures.

This guide for this skript is here: https://github.com/cortes-ciriano-lab/CancerGenomicsCourse_EMBL-EBI/tree/main

```
R --file="mut_sig_hc_filt.R" &> mut_sig_filt.Rout &
```

### Comparison to COSMIC

Following data files should be stored in your folder Project_mut_bc/data:

- COSMIC mutation data https://cancer.sanger.ac.uk/cell_lines/download:
  - CosmicCLP_MutantExport_v97_2023_01.tsv.gz
  - Cosmic_MutantCensus_v98_GRCh38.tsv.gz
- RNA edit sites, http://srv00.recas.ba.infn.it/atlas/download.html: TABLE1_hg38_chr_pos.txt
- LCRs, https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true: LCR-hs38_with_chr.bed.gz
- common SNP sites:
  - common dbSNP sites, GCF_000001405.40.common.chr.tsv
  - 1000 Genomes, 1000G_phase3_v4_20130502.sites.hg38.af01.vcf

400 verified cosmic mutations overlapping to 10 of 29 DSMZ BC cell lines were identified via the skript below. Runs a long time for loading and comparing to millions of mutation sites.

```
R --file="compare_mut2cosmic.R"
```

Get coverage/depth of theses 400 cosmic mutations for the bc cell lines via bedcov on bam files.

```
echo "Chrom Start End BT-474 CAL-120 CAL-148 CAL-51 CAL-85-1 COLO-824 DU-4475 EFM-19 EVSA-T MFM-223" > tables/cosmic_samples.csv
    samtools bedcov tables/cosmic_mut_Project_bc_mut_gatk_hc_dbsnp2.0based.bed bam/BT-474.recal.bam bam/CAL-120.recal.bam bam/CAL-148.recal.bam bam/CAL-51.recal.bam bam/CAL-85-1.recal.bam bam/COLO-824.recal.bam bam/DU-4475.recal.bam bam/EFM-19.recal.bam bam/EVSA-T.recal.bam bam/MFM-223.recal.bam >> tables/cosmic_samples.csv
```

Correct table col title via libreoffice and save into tables/cosmic_samples.xlsx

Skript for visualisation of COSMIC comparison to filtered variants:

```
R --file="compare_mut2cosmic2.R"
```

### Mutational burden

Visualise mutational burden as waterfall plot on the filtered mutations.

```
R --file="mut_vis_waterfall.R" &> mut_vis_waterfall.Rout  &
```
