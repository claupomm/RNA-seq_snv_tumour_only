###################################################
## visualisation of snp variation
## via GenVisR
## input: maf files
###################################################

library(GenVisR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readxl)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)
library(gridExtra)

set.seed(383)

expname = "_Project_bc_mut_gatk_hc_dbsnp_cosmic"
# expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp_cosmic", sep="")

# get filtered mutations with indels
# get subtypes
samples = read.csv("file2sample.csv")
samples = samples[!duplicated(samples$cell_line), ]
# Function to load vcf from haplotypecaller
Load_vcf_hc <- function(vcf_file, sample){
    VCF <- as.data.frame(fread(cmd = paste("zgrep -v '##' ", vcf_file, sep = ''), sep = '\t'))
    VCF <- dplyr::rename(VCF, CHROM = `#CHROM`)
    VCF <- VCF[VCF$CHROM %in% c(paste('chr',seq(1,22), sep = ''),'chrX','chrY'),]

    # We keep only the PASS mutations
    VCF2 <- VCF[VCF$FILTER == 'PASS',]

    # Extract tumour alternative counts
    VCF2 <- tidyr::separate(data = VCF2, col = sample, sep = ':', 
                            into = c('T_GT','T_AD','Depth','T_GQ','T_PL'))
    VCF2$REF_count <- apply(VCF2, 1, function(x) as.numeric(unlist(strsplit(x = x["T_AD"], split = ','))[1]))
    VCF2$ALT_count <- apply(VCF2, 1, function(x) as.numeric(unlist(strsplit(x = x["T_AD"], split = ','))[2]))
    VCF2$Depth <- as.numeric(VCF2$Depth) # We want this value to be considered as a number

    # Too many columns for now, let's keep only the important ones
    VCF3 <- VCF2[,c('CHROM','POS','REF','ALT','QUAL','FILTER', 'Depth','REF_count','ALT_count')]

    # Compute Variant Allele Frequency (VAF) in the tumour sample
    VCF3$VAF <- VCF3$ALT_count/VCF3$Depth
    return(VCF3)
}
# breast cancer
Mut_calls <- NULL
count = 1
for (sample in samples$sampleName) {
  bc = Load_vcf_hc(paste("snp/", sample, ".hc.pass2.lcr.dbsnp.vcf.gz", sep=""), sample)
  bc$sampleID = sample
  bc <- bc[bc$Depth >= 5,] # filter depth

  # take out mutations gnomad af>0.01 => in maf FILTER="common_variant"
  # + filter coding variants
  files_vep = paste("snp/", sample, ".hc.pass2.lcr.dbsnp.vcf.maf.gz", sep="")
  vep = read.csv(gzfile(files_vep), header=T, stringsAsFactors=F, sep = "\t", comment.char="#")
  vep = vep[grep("protein_coding", vep$BIOTYPE), c(1:16, 35:43, 47:75, 77:87, 93:114)]
  vep = vep[vep$HGVSp!="", ]
  vep = vep[vep$FILTER != "common_variant", ] # gnomad filter
  # => --af => phase 3 1000 genomes
  vep = vep[is.na(vep$AF) | vep$AF<0.01, ]
  vep$id = paste(vep$Chromosome, vep$Start_Position, sep="_")
  # account for deletion with pos-1, adaption to gatk hc counting
  tmp = which(vep$Variant_Type == "DEL")
  vep$id[tmp] = paste(vep$Chromosome[tmp], vep$Start_Position[tmp]-1, sep="_")

  bc$id = paste(bc$CHROM, bc$POS, sep="_")
  vep = vep[vep$id %in% bc$id,]

  Mut_calls[[count]] <- vep
  count = count + 1
}
# Usually faster than rbind in each iteration of the for loop
maf <- do.call("rbind", Mut_calls)
# select columns
sel = c("id", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp", "SYMBOL", "Hugo_Symbol", "Feature", "SIFT", "PolyPhen", "FILTER")
maf = maf[,sel]
# return to vep chromosome position for deletions
maf$id = paste(maf$Chromosome, maf$Start_Position, sep="_")

# take out mutations occuring in more than 20% of samples
# get table from mut_filter.R script
tmp = read.table("tables/snv_indel_maf_less20_Project_bc_mut_gatk_hc_dbsnp.tsv", sep="\t", header=T)
tmp$id = paste(tmp$Chromosome, tmp$Start_Position, sep="_")
tmp = tmp[order(tmp$id),]
maf2 = maf[order(maf$id),]
maf2 = maf2[maf2$id %in% tmp$id, colnames(maf2)!="id"]


# filter allowed mutation types
# Detected an invalid mutation type, valid values for MAF are: Nonsense_Mutation, Frame_Shift_Ins, Frame_Shift_Del, Translation_Start_Site, Splice_Site, Nonstop_Mutation, In_Frame_Ins, In_Frame_Del, Missense_Mutation, 5'Flank, 3'Flank, 5'UTR, 3'UTR, RNA, Intron, IGR, Silent, Targeted_Region, NA
# allowed = c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "Translation_Start_Site", "Splice_Site", "Nonstop_Mutation", "In_Frame_Ins", "In_Frame_Del", "Missense_Mutation", NA)
allowed = c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del", "Missense_Mutation", NA)
maf2 = maf2[maf2$Variant_Classification %in% allowed, ]

# get cosmic genes
gatk = read.table("tables/cosmic_mut_Project_bc_mut_gatk_hc_dbsnp2.tsv", header=TRUE, sep="\t") 
tab = maf2[maf2$SYMBOL %in% gatk$symbol, c(10,13,5)]
colnames(tab) = c("sample", "gene", "mutation")
myHierarchy <- data.table("mutation"=allowed, color=c("#bdbdbd", "#e31a1c", "#fb9a99", "#1f78b4", "#a6cee3", "#b2df8a"))
plotData = Waterfall(tab, mutationHierarchy = myHierarchy, geneMax = 50, plotA=NULL, plotB="frequency")

pdf(paste('plots/fig4a_waterfall',expname,'.pdf',sep=''), height=15,width=13)
drawPlot(plotData)
dev.off()

write.table(tab, paste('tables/waterfall',expname,'.csv',sep=''), row.names=F, quote=F, sep="\t")

## transform csv files to xls tables
setwd(paste(getwd(),'/tables',sep=''))
system('../csv2xlsNGS.pl')
system('rm *.csv')
setwd(gsub('/tables','',getwd()))


# distribution of mutation types
tab = data.frame(readxl::read_excel("tables/waterfall_Project_bc_mut_gatk_hc_dbsnp_cosmic.xls"))
colnames(tab) = c("Tumor_Sample_Barcode", "SYMBOL", "Variant_Classification")
tab = data.frame(table(tab$Variant_Classification))
colnames(tab) = c("Variant_Classification", "n")
pdf(paste("plots/fig4b_barplot_mut_types_waterfall", expname, ".pdf", sep = ""), width = 3, height = 5)
p1 <- ggplot(data = tab, aes(x = Variant_Classification, y = n)) +
    geom_bar(stat = "identity", position=position_dodge()) + #geom_jitter(width = 0.2) +
    theme_bw(base_size = 12, base_family = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Mutation types") +
    xlab("") + ylab("Number of mutations")
p1
dev.off()



