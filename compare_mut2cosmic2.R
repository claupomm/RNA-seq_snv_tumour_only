###################################################
## compare mutations with cosmic
###################################################
# after using bedtool for depth at mutation sites

library(psych)
library(tidyr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(gplots)
library(gridExtra)
library(cowplot)

expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp2", sep="")


# cell line samples
samples = read.csv("file2sample.csv")
samples = samples[!duplicated(samples$sampleName), ]

# get mutation table
mut2 = read.table(paste("tables/cosmic_mut", expname, ".tsv", sep = ""), header = TRUE)
mut = read.table(paste("tables/cosmic_mut_more", expname, ".tsv", sep = ""), header = TRUE)
colnames(mut2)[1] = "cell_line"
colnames(mut2)[colnames(mut2) == "gatk"] = "found"

# read depth for cosmic variant sites
depth = readxl::read_excel("tables/cosmic_samples.xlsx")
mut2$Depth = "<=5"
for (i in 1:nrow(mut2)) {
    tmp = depth[depth$Chrom %in% mut2$Chrom[i] & depth$End %in% mut2$Pos[i], match(mut2$cell_line[i], colnames(depth))]
    tmp = unique(tmp)
    if(i%%50 ==0) print(paste(date(), i))
    if (tmp>4) {
        mut2$Depth[i] = ">5"
    }
}
rm(depth)


# write table
mut2$Mutation_AA = mut$Mutation_AA[match(mut2$HGVSG,mut$HGVSG)]
mut2$Mutation_Description = mut$Mutation_Description[match(mut2$HGVSG,mut$HGVSG)]
write.table(mut2, paste("tables/cosmic_mut", expname, ".tsv", sep = ""),
    row.names = FALSE, sep = "\t")

# visualisation 
mut_long <- gather(mut2[mut2$found=="no", colnames(mut2) %in% c("cell_line", "LCR", "Freq", "Depth")], cosmic, count, LCR:Depth, factor_key = TRUE)
mut_long$cell_line = as.factor(mut_long$cell_line)
mut_long$count = factor(mut_long$count, levels=c("in", "out", ">20%", "<=20%", "<=5", ">5"))
p2 <- ggplot(data = mut_long, aes(x = cell_line, fill = count)) +
    geom_bar(stat = "count", position = "fill") + 
    theme_bw(base_size = 12, base_family = "") +
    scale_fill_brewer(palette = "Paired") +
    facet_wrap(~cosmic, nrow = 3) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.title=element_blank()) +
    ggtitle("Cosmic variants") +
    xlab("") + ylab("Portion of mutations")

# overwrite overview plot for manuscript
p1 <- ggplot(data = mut2, aes(x = cell_line, fill = found)) +
    geom_bar(stat = "count", position = position_dodge()) + 
    theme_bw(base_size = 12, base_family = "") +
    scale_fill_brewer(palette = "Paired") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Cosmic variants") +
    xlab("") + ylab("Number of mutations")


# heatmap of cosmic gene expression
tpm = read.csv("../tables/TPM_salmon_Project_bc.tsv", sep = "\t")
colnames(tpm) = gsub("\\.", "-", colnames(tpm))
genes = unique(mut2$symbol)
tab = tpm[match(genes, tpm$external_gene_name), c(9:ncol(tpm))]
rownames(tab) = genes
tab = log2(tab[!is.na(tab[,1]),]+1)
mean = apply(tab,1,mean)
tab = tab[order(mean, decreasing=F), colnames(tab) %in% unique(mut2$cell_line)]
tab_long = gather(cbind(symbol=rownames(tab),tab), cell_line, value, "BT-474":"MFM-223") 

mid = median(tab_long$value)
p3 = ggplot(tab_long, aes(cell_line, symbol)) + 
    geom_tile(aes(fill = value)) +
    theme_bw(base_size = 12, base_family = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
    labs(x = "", y = "") +
    ggtitle("Cosmic genes") +
    scale_y_discrete(limits = unique(tab_long$symbol)) +
    scale_fill_gradient2(midpoint=mid, low="#377eb8", mid="white",
                     high="#e41a1c", space ="Lab")
  
pdf(paste('plots/cosmic_heatmap_variant',expname,'.pdf',sep=''), height = 4, width=10)
plot_grid(p1, p2, p3, labels = c("A", "B", "C"), ncol=3)
dev.off()
