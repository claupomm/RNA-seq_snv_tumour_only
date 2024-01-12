###################################################
## visualise statistics
###################################################

library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cowplot)


expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp", sep="")



# load statistics for all mutation and further filter
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
all_long <- gather(all, filter, count, all:dbsnp, factor_key = TRUE)



# load mapping statistics
common_mut = fread(paste('tables/snv_indel_maf',expname,'.tsv',sep=''))
sample20 = fread(paste('tables/snv_indel_maf_less20',expname,'.tsv',sep=''))
tab = common_mut %>% count(sampleID) %>% as.data.frame()
tab = merge(tab, by="sampleID", sample20 %>% count(sampleID) %>% as.data.frame())
colnames(tab) = c("sampleID", "snv_indel", "snv_indel_20")
stats = t(read.table("tables/star_stat_summary.tsv", sep="\t", header=TRUE))
tab$mapped = stats[,"unique"]/10^6
all = cbind(all,tab[,-1])
all = all[,c(1,10, 2:9)]
all_long2 <- gather(all, filter, count, all:snv_indel_20, factor_key = TRUE)



# visualise
pdf(paste("plots/hc_stats_merged_single_violin", expname, ".pdf", sep = ""), width = 9, height = 9)
p <- ggplot(data = all_long[!all_long$filter %in% c("kG", "dbsnp", "snv_indel", "snv_indel_20") ,], aes(x = cell_line, y = count, fill = filter)) +
    geom_bar(stat = "identity", position=position_dodge()) +
    theme_bw(base_size = 12, base_family = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position="none") +
    facet_wrap(~filter, ncol = 1) +
    scale_fill_brewer(palette = "Paired") +
    ggtitle("Mutation filter") +
    xlab("") + ylab("Number of mutations")
p1 <- ggplot(data = all_long2[!all_long$filter %in% "kG",], aes(x = filter, y = count)) +
    geom_violin() + #geom_jitter(width = 0.2) +
    geom_boxplot(width=0.2, color="grey", alpha=0.2) +
    scale_y_log10() +
    theme_bw(base_size = 12, base_family = "") +
    # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Filtering steps") +
    xlab("") + ylab("Number of mutations")
# number of filtered variants correlation to star statistics
p2 <- ggplot(data = all_long2[all_long2$filter %in% c("all", "pass", "dbsnp", "snv_indel_20"),], aes(x = mapped, y = count)) +
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

