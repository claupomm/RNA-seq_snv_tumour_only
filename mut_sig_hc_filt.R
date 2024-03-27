###################################################
## mutational signatures
## restrict on SBS
###################################################

# Loading required packages
library(data.table)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(GenomicRanges)
library(MutationalPatterns)
library(plyr)
library(NMF)
library(gridExtra)
library(cowplot)


expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp", sep="")

# get subtypes
samples = read.csv("file2sample.csv")
samples = samples[!duplicated(samples$cell_line), ]



# Function to load vcf from haplotypecaller
Load_vcf_hc <- function(vcf_file, sample){
    VCF <- as.data.frame(fread(cmd = paste("zgrep -v '##' ", vcf_file, sep = ''), sep = '\t'))
    VCF <- dplyr::rename(VCF, CHROM = `#CHROM`)
    VCF <- VCF[VCF$CHROM %in% c(paste('chr',seq(1,22), sep = ''),'chrX','chrY'),]

    # keep PASS mutations
    VCF2 <- VCF[VCF$FILTER == 'PASS',]

    # Extract tumour alternative counts
    VCF2 <- tidyr::separate(data = VCF2, col = sample, sep = ':', 
                            into = c('T_GT','T_AD','T_DP','T_GQ','T_PL'))
    VCF2$T_REF_COUNT <- apply(VCF2, 1, function(x) as.numeric(unlist(strsplit(x = x["T_AD"], split = ','))[1]))
    VCF2$T_ALT_COUNT <- apply(VCF2, 1, function(x) as.numeric(unlist(strsplit(x = x["T_AD"], split = ','))[2]))
    VCF2$T_DP <- as.numeric(VCF2$T_DP) # We want this value to be considered as a number

    # Too many columns for now, let's keep only the important ones
    VCF3 <- VCF2[,c('CHROM','POS','REF','ALT','QUAL','FILTER', 'T_DP','T_REF_COUNT','T_ALT_COUNT')]

    # extract snv only, filter out indels
    VCF3 = VCF3[apply(VCF3[, 3:4], 1, function(x) nchar(x[1])==1 & nchar(x[2])==1),]

    # Compute Variant Allele Frequency (VAF) in the tumour sample
    VCF3$T_VAF <- VCF3$T_ALT_COUNT/VCF3$T_DP
    return(VCF3)
}


##############################
## add own data, bc
##############################
# breast cancer
Mut_calls <- NULL
count = 1
for (sample in samples$sampleName) {
  bc = Load_vcf_hc(paste("snp/", sample, ".hc.pass2.lcr.dbsnp.vcf.gz", sep=""), sample)
  bc$sampleID = sample
  bc <- bc[bc$T_DP >= 5,] # filter depth

  # take out mutations gnomad af>0.01 => in maf FILTER="common_variant"
  # + filter coding variants
  files_vep = paste("snp/", sample, ".hc.pass2.lcr.dbsnp.vcf.maf.gz", sep="")
  vep = read.csv(gzfile(files_vep), header=T, stringsAsFactors=F, sep = "\t", comment.char="#")
  vep = vep[grep("protein_coding", vep$BIOTYPE), c(1:16, 35:43, 47:75, 77:87, 93:114)]
  vep = vep[vep$HGVSp!="", ]
  vep = vep[vep$FILTER != "common_variant", ] # gnomad filter
  # => --af => phase 3 1000 genomes
  vep = vep[is.na(vep$AF) | vep$AF<0.01, ]

  vep$id = paste(vep$Chromosome, vep$Start_Position, vep$Reference_Allele, vep$Tumor_Seq_Allele2, sep="_")
  bc$id = paste(bc$CHROM, bc$POS, bc$REF, bc$ALT, sep="_")
  bc = bc[bc$id %in% vep$id,]

  Mut_calls[[count]] <- bc
  count = count + 1
}
# Usually faster than rbind in each iteration of the for loop
Mut_calls <- do.call("rbind", Mut_calls)


# filter mutations occuring in more than 20% of samples
Mut_calls$chr_start = paste(Mut_calls$CHROM, Mut_calls$POS)
Mut_calls = Mut_calls[order(Mut_calls$chr_start),]
write.table(Mut_calls, paste('tables/snv_maf',expname,'.csv',sep=''), row.names=F, quote=F, sep="\t")
tmp = as.data.frame(table(Mut_calls$chr_start))
tmp = as.vector(tmp$Var1[tmp$Freq >= ceiling(nrow(samples)/5)])
Mut_calls = Mut_calls[!Mut_calls$chr_start %in% tmp, colnames(Mut_calls)!="chr_start"]
Mut_calls = Mut_calls[order(Mut_calls$sampleID),]
write.table(Mut_calls, paste('tables/snv_maf_less20',expname,'.csv',sep=''), row.names=F, quote=F, sep="\t")
# Mut_calls = readxl::read_excel("tables/snv_maf_less20_Project_bc_mut_gatk_hc_dbsnp.xls")



## transform csv files to xls tables
setwd(paste(getwd(),'/tables',sep=''))
system('../csv2xlsNGS.pl')
system('rm *.csv')
setwd(gsub('/tables','',getwd()))



# Mutation characteristics
# Transform data.frame to genomic ranges object
# Create list of granges
GRmaker <- function(CALLS,GROUP){
  TEMP <- CALLS[CALLS$GROUP == GROUP,]
  GR <- with(TEMP, GRanges(CHROM, ranges=IRanges(POS, POS), REF = REF,ALT = ALT))
  seqlevelsStyle(GR) <- "UCSC"
  GenomeInfoDb::genome(GR) <- 'hg38'
  return (GR)
}


# Once defined this function, we can execute it using our mutations data frame:
# This column will specify how to group the mutations
# For now, we use the sampleID, but it can be modified
Mut_calls$GROUP <- Mut_calls$sampleID

# Create list of granges
GROUPS <- unique(Mut_calls$GROUP)
G_lists <-  lapply(GROUPS, function(x) GRmaker(Mut_calls, x))
names(G_lists) <- GROUPS


# Getting known signatures
# breast cancer mut sig, Alexandrov 2020
signatures = signatures[,colnames(signatures) %in% paste("SBS", c(1,2,3,5,8,9,13,"17a", "17b",18,37,40,41), sep="")]


# The function fit_to_signatures_strict calculates the best possible linear combination of mutational signatures, adhering to non-negative least-squares constraints, in order to reconstruct the mutation matrix with the highest degree of accuracy.
# max_delta parameter: Decreasing this number will make the refitting less strict, while increasing it will make the refitting more strict. default 0.004
strict_refit <- fit_to_signatures_strict(mut_mat, as.matrix(signatures), max_delta = 0.001, method = "best_subset")

# Getting the contribution of each signature
fit_res_strict <- strict_refit$fit_res


pdf(paste("plots/mut_sig", expname, ".pdf", sep = ""), width = 9, height = 9)
type_occurrences <- mut_type_occurrences(G_lists, ref_genome)
p1 <- plot_spectrum(type_occurrences, indv_points = TRUE)
p2 <- plot_96_profile(mut_mat[,c(1,4, 6,9,15, 25,28,29)]) # cell lines with typical mut sig
p3 <- plot_contribution(fit_res_strict$contribution,
        coord_flip = TRUE,
        mode = "relative"
)
p13 = plot_grid(p1, p3, labels = c("A", "C"), ncol=1, rel_heights = c(1, 2))
plot_grid(p13, p2, ncol=2, labels = c("", "B"), rel_widths = c(1, 1))
dev.off()




sessionInfo()

