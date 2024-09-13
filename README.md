# Variant calling pipeline on RNA-seq data for tumour-only breast cancer cell lines

Typically, variant calling on tumour material is performed on whole-exome or whole-genome DNA sequencing. This workflow, however, has been used for analysing RNA-seq data of 29 tumour breast cancer cell lines without matched-normal sample pairs for identifying single nuleotide variants (SNVs) and small insertions / deletions (InDels). Since capturing germline variants produces a massive amount of mutations per sample, the filtering process is of special importance in order to retrieve an essence of potentially intriguing variants. \
This proposed pipeline could serve as a guideline how to proceed with RNA-seq tumour-only samples and probably needs adjustments for each project depending on the tumour type and the user's need and scope. \
Please note, that third-party data download links might change over time and so table columns or data formats, which also would require script adaption.



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
5. Variant validation
	1. Statistics
	2. Mutational signatures
	3. Comparison to COSMIC
	4. Mutational burden
	5. Specificity and sensitivity



## Intial steps


### Tools

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


### Data

Following data need to be downloaded to your project folder to the subfolder data, some websites need registration:

- for GATK:
  
  - sources: 
    
    - https://gatk.broadinstitute.org/hc/en-us/articles/360035890811
    - or look here: https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&supportedpurview=project&prefix=&forceOnObjectsSortingFiltering=false
  
  - files 
    
    - Homo_sapiens_assembly38.fasta
    - 1000G_phase1.snps.high_confidence.hg38.vcf
    - 1000G_phase3_v4_20130502.sites.hg38.vcf
    - Mills_and_1000G_gold_standard.indels.hg38.vcf

- RNA edit sites, http://srv00.recas.ba.infn.it/atlas/download.html: TABLE1_hg38_chr_pos.txt

- LCRs, https://github.com/lh3/varcmp/blob/master/scripts/LCR-hs38.bed.gz?raw=true: LCR-hs38_with_chr.bed.gz (curl download see below)

- COSMIC mutation data https://cancer.sanger.ac.uk/cosmic/download/cell-lines-project#cell-lines-project:
  
  - CosmicCLP_MutantExport_v97_2023_01.tsv.gz
  - Cosmic_MutantCensus_v98_GRCh38.tsv.gz

- Common SNP sites
  
  - dbSNP, GCF_000001405.40.common.chr.tsv, https://ftp.ncbi.nih.gov/snp/latest_release/VCF/
  - 1000 Genomes, 1000G_phase3_v4_20130502.sites.hg38.af01.vcf
  - gnomAD: contained in VEP data sources
  

- COSMIC mutational signature data:
  
	+ https://cancer.sanger.ac.uk/signatures/downloads/
	+ COSMIC_catalogue-signatures_SBS96_v3.3.zip
	+ COSMIC_v3.3.1_SBS_GRCh38.txt
  

General pipeline scheme:
![Overview flowchart](./plots/Figure_S2_flowchart.pdf "Overview flowchart")


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

These should be downloaded to your directory /path/to/fastq_files. Afterwards these raw fastq files are linked and prepared for the pipeline:

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


### Required R libraries

For analysis and visualisation following R packages are required:
- BSgenome, BSgenome.Hsapiens.UCSC.hg38, cowplot, data.table, dplyr, ggplot2, ggpubr, gplots, gridExtra, gtools, GenomicRanges, GenVisR, MutationalPatterns, NMF, plyr, psych, RColorBrewer, tidyr, readxl

You can install the package within R (version 4.4.1) like so:
```
BiocManager::install(c("BSgenome", "BSgenome.Hsapiens.UCSC.hg38", "cowplot", "data.table", "dplyr", "ggplot2", "ggpubr", "gplots", "gridExtra", "gtools", "GenomicRanges", "GenVisR", "MutationalPatterns", "NMF", "plyr", "psych", "RColorBrewer", "tidyr", "readxl"))
```



## Preprocessing, mapping, and variant calling

The preprocessing steps are best to be done on a high computing cluster. Here, we do this via the job scheduler SLURM.

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

Create an overview on the sequencing depth for each sample:

```
R --file="stats_seq.R"
```



## Filtering steps


### Filter RNA-edit sites

Some basic knowledge about RNA-edit sites is explained here: https://academic.oup.com/nar/article/49/D1/D1012/5940507?login=true

```
rna_edit=data/TABLE1_hg38_chr_pos.txt
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


### Filter common SNPs from called SNVs

Prepare the dbSNP filter together with the downloaded dbsnp file at https://ftp.ncbi.nih.gov/snp/latest_release/VCF/ stored in $DIR/data

```
cd $DIR/data

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
zgrep -c "^chr" *.hc.pass2.lcr.1kG.vcf.gz > stats_hc_1kG.txt # q20 + depth + rna_edit + lcr + 1000Gzgrep -c "^chr" *.hc.pass2.lcr.dbsnp.vcf.gz > stats_hc_dbsnp.txt # q20 + depth + rna_edit + lcr + dbsnp
```

Visualise variant number after each step of filtering:

```
R --file="stats.R"
```


### Comparison to COSMIC

Following data files should be stored in your folder Project_mut_bc/data:

- COSMIC mutation data https://cancer.sanger.ac.uk/cosmic/download/cell-lines-project#cell-lines-project (needs login):
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

Following cosmic_samples.xlsx was already provided, hence can be skipped and you may proceed to compare_mut2cosmic2.R, but if wanted, you may create the coverage table yourself:

Get coverage/depth of theses 400 cosmic mutations for the bc cell lines via bedcov on bam files.

```
echo -e "Chrom\tStart\tEnd\tBT-474\tCAL-120\tCAL-148\tCAL-51\tCAL-85-1\tCOLO-824\tDU-4475\tEFM-19\tEVSA-T\tMFM-223" > tables/cosmic_samples.csv
samtools bedcov tables/cosmic_mut_Project_bc_mut_gatk_hc_dbsnp2.0based.bed bam/BT-474.recal.bam bam/CAL-120.recal.bam bam/CAL-148.recal.bam bam/CAL-51.recal.bam bam/CAL-85-1.recal.bam bam/COLO-824.recal.bam bam/DU-4475.recal.bam bam/EFM-19.recal.bam bam/EVSA-T.recal.bam bam/MFM-223.recal.bam >> tables/cosmic_samples.csv
```

Open table and save into tables/cosmic_samples.xlsx 

Skript for visualisation of COSMIC comparison to filtered variants:

```
R --file="compare_mut2cosmic2.R"
```


### Mutational signatures

Filter mutations by low complexity regions LCR, gnomad af > 0.01, extract single base substitutions (SBS). The R-package MutationalPatterns is used for calculating mutational signatures.

A guide for this skript is here: https://github.com/cortes-ciriano-lab/CancerGenomicsCourse_EMBL-EBI/tree/main

```
R --file="mut_sig_hc_filt.R" &> mut_sig_filt.Rout &
```


### Mutational burden

Visualise mutational burden as waterfall plot on the filtered mutations.

```
R --file="mut_vis_waterfall.R" &> mut_vis_waterfall.Rout  &
```


### Specificity and sensitivity

In order to compare the improvement of the successive filtering steps of the pipeline, specificity and sensitivity are calculated based on the verified COSMIC 400 variants set of the 10 overlapping breast cancer cell lines of DSMZ. 
Counting of the variant intersect is carried out via vcftools (https://vcftools.github.io/). 

#### Prepare sample vcf
As vcftools needs the vcf format for comparison, the COSMIC variants are extracted from the corresponding unfiltered called variant vcf for each cell line and stored as vcf for vcftools.
```
cd $DIR
R --file="spec_sens.R"
```

#### Identify intersection and unique variants
For each sample and filter step the overlapping and unique variants to the defined COSMIC variant set are detected via vcftools. 

```
cd $DIR
start="cosmic_mut_more_chr_pos_"
end="_Project_bc_mut_gatk_hc_dbsnp2.vcf"
for i in $(cut -d, -f1 file2sample_cosmic.csv| tail -n+2 - ); do
date
echo $i
chr=$(cut -f1 spec_sens/${start}$i${end} | grep "^chr" | uniq | tr '\n' ',')
chr="${chr/%,/}"
chr="${chr/chr/--chr chr}"
chr="${chr//,/ --chr }"
vcftools --diff-site --out "spec_sens/$i.cosmic400.hc" $chr --gzvcf snp/$i.hc.pass.vcf.gz --diff spec_sens/${start}$i${end}
vcftools --diff-site --out "spec_sens/$i.cosmic400.hc.pass" $chr --gzvcf snp/$i.hc.pass.vcf.gz --gzdiff spec_sens/${start}$i${end}
vcftools --diff-site --out "spec_sens/$i.cosmic400.hc.pass2" $chr --gzvcf snp/$i.hc.pass2.vcf.gz --gzdiff spec_sens/${start}$i${end}
vcftools --diff-site --out "spec_sens/$i.cosmic400.hc.pass2.lcr" $chr --gzvcf snp/$i.hc.pass2.lcr.vcf.gz --gzdiff spec_sens/${start}$i${end}
vcftools --diff-site --out "spec_sens/$i.cosmic400.hc.pass2.lcr.dbsnp" $chr --gzvcf snp/$i.hc.pass2.lcr.dbsnp.vcf.gz --gzdiff spec_sens/${start}$i${end}
done
```

#### Summarise
Count all unique/overlapping variants for each filter step and sample to one table.
```
cd $DIR/spec_sens
echo "sample,TP_hc,TP_pass,TP_pass2,TP_lcr,TP_dbsnp,FP_hc,FP_pass,FP_pass2,FP_lcr,FP_dbsnp,FN_hc,FN_pass,FN_pass2,FN_lcr,FN_dbsnp" > overlap_cosmic400_samples.csv
end=diff.sites_in_files
for j in $(cut -d, -f1 ../file2sample_cosmic.csv| tail -n+2 - ); do
echo $j
i=$j.cosmic400
echo $j","$(cut -f4 $i.hc.$end | grep -c B)","$(cut -f4 $i.hc.pass.$end | grep -c B)","$(cut -f4 $i.hc.pass2.$end | grep -c B)","$(cut -f4 $i.hc.pass2.lcr.$end | grep -c B)","$(cut -f4 $i.hc.pass2.lcr.dbsnp.$end | grep -c B)","$(cut -f4 $i.hc.$end | grep -c 1)","$(cut -f4 $i.hc.pass.$end | grep -c 1)","$(cut -f4 $i.hc.pass2.$end | grep -c 1)","$(cut -f4 $i.hc.pass2.lcr.$end | grep -c 1)","$(cut -f4 $i.hc.pass2.lcr.dbsnp.$end | grep -c 1),"$(cut -f4 $i.hc.$end | grep -c 2)","$(cut -f4 $i.hc.pass.$end | grep -c 2)","$(cut -f4 $i.hc.pass2.$end | grep -c 2)","$(cut -f4 $i.hc.pass2.lcr.$end | grep -c 2)","$(cut -f4 $i.hc.pass2.lcr.dbsnp.$end | grep -c 2)" >> overlap_cosmic400_samples.csv
done
```

#### Calculate and visualise
As basis for calculating sensitivity and specificity the numbers of all positive variants of a sample are derived from the extracted 400 COSMIC variants and the numbers of all negative variants are taken from the unfiltered called variants of a sample without the positive variants. Furtheron, true positive variants correspond to the intersect to the COSMIC variant set and false positive are those, which are unique and not overlapping to the COSMIC variant set for each sample and filter step, and false negatives those, which are unique for COSMIC variant set.
```
cd $DIR
R --file="spec_sens2.R"
```