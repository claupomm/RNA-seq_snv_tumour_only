###################################################################
# Calculate and visualise sensitivity, specificity, and precision:
# As basis for calculating sensitivity, specificity, and precision 
# the numbers of all positive variants of a sample are derived 
# from the extracted 400 COSMIC variants and 
# the numbers of all negative variants are taken 
# from the unfiltered called variants of a sample without the positive variants. 
# Furtheron, true positive variants correspond to 
# the intersect to the COSMIC variant set and false positive are those, 
# which are unique and not overlapping to the COSMIC variant set 
# for each sample and filter step, 
# and false negatives those, which are unique for COSMIC variant set.
###################################################################

# filter steps:
# 1. identify SNV via haplotypecaller => .hc.vcf.gz
# 2. <5 depth => .hc.pass.vcf.gz
# 3. RNA edit sites => .hc.pass2.vcf.gz
# 4. lcr sites => .hc.pass2.lcr.vcf.gz
# 5. dbsnp => .hc.pass2.lcr.dbsnp.vcf.gz
# 6. not protein coding + 1000genomes + gnomAD => R, snv_indel
# 7. >=20% => R, snv_indel_20

library(psych)
library(tidyr)
library(ggplot2)
library(data.table)
library(cowplot)
library(gtools)

expname = "_Project_bc_mut_gatk_hc_dbsnp2"
# expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp2", sep="")


# cell line samples
samples = read.csv("file2sample.csv")
samples = samples[!duplicated(samples$sampleName), ]


# get results from filtering dp5/gnomad/1000G/dbsnp/protein-coding
ref = read.table(paste("tables/cosmic_mut", expname, ".tsv", sep = ""),
    header=T, sep = "\t")
start="cosmic_mut_more_chr_pos_"
end="_Project_bc_mut_gatk_hc_dbsnp2.vcf"

spec = read.table("spec_sens/overlap_cosmic400_samples.csv", header=T, sep = ",", row.names=1)
colnames(spec) = gsub("_pass2", "_edit", colnames(spec))
spec$P_cosmic_all = 0
spec$N = spec$FP_hc

for (i in unique(ref$cell_line)) {
    tmp = ref[ref$cell_line==i,]
    tmp = tmp[order(tmp$Pos),]
    tmp = tmp[mixedorder(tmp$Chrom),]
    # print(paste("duplicated:", sum(duplicated(paste(tmp$Chrom, tmp$Pos)))))
    spec[i, "P_cosmic_all"] = nrow(tmp)
}


# calculate specificity and sensitivity for the different mutation filtering steps
# sensitivity: TP/(TP+FN)
# specificity: TN/(TN+FP)
# precision: TP/(TP+FP)
tab = data.frame(cell_line = rownames(spec))
i = "hc"
tmp = cbind(spec[,paste("TP", i, sep="_")], 
            spec[,paste("FP", i, sep="_")],
            spec[,paste("FN", i, sep="_")],
            spec[,"P_cosmic_all"],
            spec[,"N"])
tmp = data.frame(tmp)
tmp = cbind(rownames(spec), filter=i, tmp)
colnames(tmp) = c("cell_line", "filter", "TP", "FP", "FN", "P", "N")
tmp$specificity = (tmp$N-tmp$FP)/tmp$N
tmp$sensitivity = tmp$TP/tmp$P
tmp$precision = tmp$TP/(tmp$TP+tmp$FP)
tab = tmp
for (i in c("pass", "edit", "lcr", "dbsnp")){
    tmp = cbind(spec[,paste("TP", i, sep="_")], 
            spec[,paste("FP", i, sep="_")],
            spec[,paste("FN", i, sep="_")],
            spec[,"P_cosmic_all"],
            spec[,"N"])
    tmp = data.frame(tmp)
    tmp = cbind(rownames(spec), filter=i, tmp)
    colnames(tmp) = c("cell_line", "filter", "TP", "FP", "FN", "P", "N")
    tmp$sensitivity = tmp$TP/tmp$P
    tmp$specificity = (tmp$N-tmp$FP)/tmp$N
    tmp$precision = tmp$TP/(tmp$TP+tmp$FP)
    tab = rbind(tab, tmp)
}
tab$filter = factor(x = tab$filter, levels = unique(tab$filter))



pdf(paste("plots/fig_s1_spec_sens", expname, ".pdf", sep = ""), width = 9, height = 3)
p <- ggplot(data = tab, aes(x = filter, y = sensitivity)) +
    geom_violin(width=0.2) +
    geom_point(color="grey") +
    ylim(0, 1) +
    theme_bw(base_size = 12, base_family = "") +
    ggtitle("Sensitivity after filtering") +
    xlab("") + ylab("Sensitivity")
p1 <- ggplot(data = tab, aes(x = filter, y = specificity)) +
    geom_violin(outlier.shape = NA) +
    geom_boxplot(width=0.2, color="grey", alpha=0.2) +
    ylim(0, 1) +
    theme_bw(base_size = 12, base_family = "") +
    ggtitle("Specificity after filtering") +
    xlab("") + ylab("Specificity")
p2 <- ggplot(data = tab, aes(x = filter, y = precision)) +
    geom_violin(outlier.shape = NA) +
    geom_boxplot(width=0.2, color="grey", alpha=0.2) +
    # ylim(0, 1) +
    theme_bw(base_size = 12, base_family = "") +
    ggtitle("Precision after filtering") +
    xlab("") + ylab("Precision")
plot_grid(p, p1, p2, ncol=3, labels = c("a", "b", "c"))
dev.off()