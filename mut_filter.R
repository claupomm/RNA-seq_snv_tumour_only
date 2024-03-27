###################################################
## get mutation numbers for each sample
## after filtering common mutants
## gnomad, dbsnp, 1000genomes
## and after >20% sample filter
###################################################

library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp", sep="")


# get filtered mutations with indels
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

    # Extract maf file
    VCF2 <- tidyr::separate(data = VCF2, col = sample, sep = ':', 
                            into = c('T_GT','T_AD','Depth','T_GQ','T_PL'))
    VCF2$REF_count <- apply(VCF2, 1, function(x) as.numeric(unlist(strsplit(x = x["T_AD"], split = ','))[1]))
    VCF2$ALT_count <- apply(VCF2, 1, function(x) as.numeric(unlist(strsplit(x = x["T_AD"], split = ','))[2]))
    VCF2$Depth <- as.numeric(VCF2$Depth) # We want this value to be considered as a number

    # Keep interesting columns
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
  # account for deletion with pos-1
  tmp = which(vep$Variant_Type == "DEL")
  vep$id[tmp] = paste(vep$Chromosome[tmp], vep$Start_Position[tmp]-1, sep="_")

  bc$id = paste(bc$CHROM, bc$POS, sep="_")
  bc = bc[bc$id %in% vep$id,]
  bc = merge(bc, vep[,c("id", "Chromosome", "Start_Position", "End_Position", "Strand", "HGVSc", "HGVSp", "HGVSp_Short", "Variant_Classification", "Consequence", "Existing_variation", "Hugo_Symbol", "Gene", "Entrez_Gene_Id", "SIFT", "PolyPhen", "IMPACT")], by="id")

  Mut_calls[[count]] <- bc
  count = count + 1
}
# Usually faster than rbind in each iteration of the for loop
Mut_calls <- do.call("rbind", Mut_calls)
# rearrange column order
Mut_calls = Mut_calls[,c(1,12:16, 4:5, 8:11, 17:ncol(Mut_calls))]

# filter mutations occuring in more than 20% of samples
tmp = as.data.frame(table(Mut_calls$id))
tmp = as.vector(tmp$Var1[tmp$Freq >= ceiling(nrow(samples)/5)])
Mut_calls$Freq = "<=20%"
Mut_calls$Freq[Mut_calls$id %in% tmp] = ">20%"
write.table(Mut_calls[,-1], paste('tables/snv_indel_maf',expname,'.tsv',sep=''), row.names=F, quote=F, sep="\t")
write.table(Mut_calls[!Mut_calls$id %in% tmp,-1], paste('tables/snv_indel_maf_less20',expname,'.tsv',sep=''), row.names=F, quote=F, sep="\t")

# merge summary in one table
common_mut = data.table::fread(paste('tables/snv_indel_maf',expname,'.tsv',sep=''))
sample20 = data.table::fread(paste('tables/snv_indel_maf_less20',expname,'.tsv',sep=''))
tab = common_mut %>% count(sampleID) %>% as.data.frame()
tab = merge(tab, by="sampleID", sample20 %>% count(sampleID) %>% as.data.frame())
colnames(tab) = c("sampleID", "snv_indel", "snv_indel_20")
write.table(tab, paste('tables/snv_indel_summary',expname,'.tsv',sep=''), row.names=F, quote=F, sep="\t")

# violin plot
tab = read.table(paste('tables/snv_indel_summary',expname,'.tsv',sep=''), header=TRUE)
all_long <- gather(tab, filter, count, snv_indel:snv_indel_20, factor_key = TRUE)
pdf(paste("plots/hc_stats_violin", expname, ".pdf", sep = ""), width = 3, height = 3)
p <- ggplot(data = all_long[grep("_p$", all_long$filter, invert=TRUE),], aes(x = filter, y = count)) +
    geom_violin(outlier.shape = NA) + #geom_jitter(width = 0.2) +
    geom_boxplot(width=0.2, color="grey", alpha=0.2) +
    scale_y_log10(limits=c(10^2,10^4)) +
    theme_bw(base_size = 12, base_family = "") +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    # facet_wrap(~filter, nrow = 3) + 
    # scale_fill_brewer(palette = "Paired") +
    ggtitle("In coding regions") +
    xlab("") + ylab("Number of mutations")
p
dev.off()



# merge with number for all mutation and further filter
all = read.csv("snp/stats_hc_all.txt", sep=":", header=FALSE)
colnames(all) = c("cell_line", "all")
all$cell_line = gsub("\\.hc\\.vcf\\.gz", "", all$cell_line)
pass = read.csv("snp/stats_hc_pass.txt", sep=":", header=FALSE)
pass2 = read.csv("snp/stats_hc_pass2.txt", sep=":", header=FALSE)
lcr = read.csv("snp/stats_hc_lcr.txt", sep=":", header=FALSE)
kG = read.csv("snp/stats_hc_1kG.txt", sep=":", header=FALSE)
dbsnp = read.csv("snp/stats_hc_dbsnp.txt", sep=":", header=FALSE)
all$pass = pass$V2
all$edit = pass2$V2
all$lcr = lcr$V2
all$kG = kG$V2
all$dbsnp = dbsnp$V2
all_long2 <- gather(all, filter, count, all:dbsnp, factor_key = TRUE)

library(gridExtra)
pdf(paste("plots/hc_stats_violin2", expname, ".pdf", sep = ""), width = 3, height = 6)
p1 <- ggplot(data = all_long2, aes(x = filter, y = count)) +
    geom_violin(outlier.shape = NA) + #geom_jitter(width = 0.2) +
    geom_boxplot(width=0.2, color="grey", alpha=0.2) +
    scale_y_log10(limits=c(10^4,10^6)) +
    theme_bw(base_size = 12, base_family = "") +
    ggtitle("Filtering") +
    xlab("") + ylab("Number of mutations")
p2 <- ggplot(data = all_long, aes(x = filter, y = count)) +
    geom_violin(outlier.shape = NA) + #geom_jitter(width = 0.2) +
    geom_boxplot(width=0.2, color="grey", alpha=0.2) +
    scale_y_log10(limits=c(10^2,10^4)) +
    theme_bw(base_size = 12, base_family = "") +
    ggtitle("") +
    xlab("") + ylab("Number of mutations")
grid.arrange(p1, p2)
dev.off()

pdf(paste("plots/hc_stats", expname, ".pdf", sep = ""), width = 7, height = 7)
p <- ggplot(data = all_long2[all_long2$filter %in% c("all", "pass", "edit", "lcr") ,], aes(x = cell_line, y = count, fill = filter)) +
    geom_bar(stat = "identity", position=position_dodge()) + 
    theme_bw(base_size = 12, base_family = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    facet_wrap(~filter, nrow = 4) + 
    scale_fill_brewer(palette = "Paired") +
    ggtitle("Mutation filter") +
    xlab("") + ylab("Number of mutations")
p
dev.off()


# number of filtered variants correlation to star statistics
library(ggpubr)
common_mut = data.table::fread(paste('tables/snv_indel_maf',expname,'.tsv',sep=''))
sample20 = data.table::fread(paste('tables/snv_indel_maf_less20',expname,'.tsv',sep=''))
tab = common_mut %>% count(sampleID) %>% as.data.frame()
tab = merge(tab, by="sampleID", sample20 %>% count(sampleID) %>% as.data.frame())
colnames(tab) = c("sampleID", "snv_indel", "snv_indel_20")
stats = t(read.table("tables/star_stat_summary.tsv", sep="\t", header=TRUE))
tab$mapped = stats[,"unique"]/10^6
all = cbind(all,tab[,-1])
all = all[,c(1,10, 2:9)]
all_long <- gather(all, filter, count, all:snv_indel_20, factor_key = TRUE)
pdf(paste("plots/hc_stats_library_size_vs_variant_number", expname, ".pdf", sep = ""), width = 7, height = 7)
p <- ggplot(data = all_long, aes(x = mapped, y = count)) +
    ggtitle("Mapped reads vs variant number") +
    xlab("Mapped reads in Mio") + ylab("Number of mutations") +
    theme_bw(base_size = 12, base_family = "") +
    facet_wrap( ~ filter, ncol=3, scales="free_y") +
    stat_smooth(method = "lm", color="black", formula = y ~ x) + # linear regression line, by default includes 95% confidence region
    geom_point() + 
    stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
            r.accuracy = 0.01)
            # label.x = 30, label.y = 300, size = 4)
    # stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
    #                         label.x = 30, label.y = 400, size = 4)
p
dev.off()
pdf(paste("plots/hc_stats_library_size_vs_variant_number2", expname, ".pdf", sep = ""), width = 3, height = 10)
p <- ggplot(data = all_long[all_long$filter %in% c("all", "pass", "dbsnp", "snv_indel_20"),], aes(x = mapped, y = count)) +
    ggtitle("Mapped reads vs variant number") +
    xlab("Mapped reads in Mio") + ylab("Number of mutations") +
    theme_bw(base_size = 12, base_family = "") +
    facet_wrap( ~ filter, ncol=1, scales="free_y") +
    stat_smooth(method = "lm", color="black", formula = y ~ x) + # linear regression line, by default includes 95% confidence region
    geom_point() + 
    stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
            r.accuracy = 0.01)
            # label.x = 30, label.y = 300, size = 4)
    # stat_regline_equation(aes(label = ..eq.label..), # adds equation to linear regression
    #                         label.x = 30, label.y = 400, size = 4)
p
dev.off()

library(gridExtra)
library(cowplot)
pdf(paste("plots/hc_stats_merged_single_violin", expname, ".pdf", sep = ""), width = 9, height = 9)
p <- ggplot(data = all_long2[!all_long2$filter %in% c("kG", "dbsnp", "snv_indel", "snv_indel_20") ,], aes(x = cell_line, y = count, fill = filter)) +
    geom_bar(stat = "identity", position=position_dodge()) + 
    theme_bw(base_size = 12, base_family = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position="none") +
    facet_wrap(~filter, ncol = 1) + 
    scale_fill_brewer(palette = "Paired") +
    ggtitle("Mutation filter") +
    xlab("") + ylab("Number of mutations")
p1 <- ggplot(data = all_long[!all_long2$filter %in% "kG",], aes(x = filter, y = count)) +
    geom_violin() + #geom_jitter(width = 0.2) +
    geom_boxplot(width=0.2, color="grey", alpha=0.2) +
    scale_y_log10() +
    theme_bw(base_size = 12, base_family = "") +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Filtering steps") +
    xlab("") + ylab("Number of mutations")
p2 <- ggplot(data = all_long[all_long$filter %in% c("all", "pass", "dbsnp", "snv_indel_20"),], aes(x = mapped, y = count)) +
    ggtitle("Reads vs variant") +
    xlab("Mapped reads in Mio") + ylab("Number of mutations") +
    theme_bw(base_size = 12, base_family = "") +
    facet_wrap( ~ filter, ncol=1, scales="free_y") +
    stat_smooth(method = "lm", color="black", formula = y ~ x) + # linear regression line, by default includes 95% confidence region
    geom_point() + 
    stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
            r.accuracy = 0.01)
pp1 = plot_grid(p, p1, labels = c("A", "B"), ncol=1, rel_heights = c(2, 1))
plot_grid(pp1, p2, ncol=2, labels = c("", "C"), rel_widths = c(5, 2))
dev.off()

