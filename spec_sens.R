###################################################################
# Create vcf files for each sample for vcftools:
# As vcftools needs the vcf format for comparison, 
# the COSMIC variants are extracted from 
# the corresponding unfiltered called variant vcf for each cell line 
# and stored as vcf for vcftools.
###################################################################

# calculate specificity and sensitivity for the different mtuation filtering steps
# sensitivity: TP/(TP+FN)
# specificity: TN/(TN+FP)

# filter steps:
# 1. identify SNV via haplotypecaller => .hc.vcf.gz
# 2. >=5 depth => .hc.pass.vcf.gz
# 3. RNA edit sites => .hc.pass2.vcf.gz
# 4. lcr sites => .hc.pass2.lcr.vcf.gz
# 5. dbsnp => .hc.pass2.lcr.dbsnp.vcf.gz
# 6. protein coding + 1000genomes + gnomAD => R, snv_indel
# 7. <20% => R, snv_indel_20



library(psych)
library(tidyr)
library(ggplot2)
library(data.table)

expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp2", sep="")


# cell line samples
samples = read.csv("file2sample.csv")
samples = samples[!duplicated(samples$sampleName), ]


# get results from filtering dp5/gnomad/1000G/dbsnp/protein-coding
ref = read.table(paste("tables/cosmic_mut_more", expname, ".tsv", sep = ""),
    header=T, sep = "\t")
start="cosmic_mut_more_chr_pos_"
end="_Project_Claudia231017_bc_mut_gatk_hc_dbsnp2.vcf"

for (i in unique(ref$Sample_name)) {
    tmp = ref[ref$Sample_name==i,]
    tmp = tmp[order(tmp$Pos),]
    tmp$Chrom = as.numeric(gsub("chr", "", tmp$Chrom))
    tmp = tmp[order(tmp$Chrom),]
    tmp$Chrom[is.na(tmp$Chrom)] = "X"
    tmp$Chrom = paste("chr",tmp$Chrom, sep="")
    colnames(tmp)[colnames(tmp)=="Chrom"] = "#Chrom"
    colnames(tmp) = toupper(colnames(tmp))
    spec = paste("spec_sens/",start, i, end, sep="")
    vcf = paste("snp/",i,".hc.pass.vcf.gz", sep="")
    system(paste("zgrep ^# ", vcf, " > ", spec, sep=""))
    for (j in 1:nrow(tmp)) {
        # search for POS without last digit for fuzzy search needed for deletion/insertion
        tmp2 = system(paste('zgrep ', tmp$'#CHROM'[j],' ', vcf,' | grep -P "\t', substr(tmp$POS[j], 0, nchar(tmp$POS[j])-1),'"', sep=""), intern=T)
        if (length(tmp2)!=0) {
            system(paste('echo "', tmp2, '" >> ', spec, sep=""))
        }
    }
}

