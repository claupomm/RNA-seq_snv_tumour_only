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

expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp", sep="")

# get subtypes
samples = read.csv("file2sample.csv")
samples = samples[!duplicated(samples$cell_line), ]

# get sample mafs, gatk results
files = list.files("snp", pattern=".hc.pass2.lcr.dbsnp.vcf.maf.gz$", full.names=T,include.dirs=T)

for (i in 1:length(files)) {
	sample = gsub(".hc.pass2.lcr.dbsnp.vcf.maf.gz", "", gsub("snp/", "", files[i]))
	if(i==1) {
		maf = read.csv(gzfile(files[i]), sep="\t", skip=1)
		maf$Tumor_Sample_Barcode = sample
	} else {
		tmp = read.csv(gzfile(files[i]), sep="\t", skip=1)
		tmp$Tumor_Sample_Barcode = sample
		maf = rbind(maf,tmp)
	}
}


# Focus on specific mutation types
allowed = c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "Translation_Start_Site", "Splice_Site", "Nonstop_Mutation", "In_Frame_Ins", "In_Frame_Del", "Missense_Mutation", NA)

# filter common mutations
# => --af => phase 3 1000 genomes
maf2 = maf[maf$Variant_Classification %in% allowed, ]
maf2 = maf2[is.na(maf2$AF) | maf2$AF<0.01, ]

# get cosmic genes
gatk = read.table("tables/cosmic_mut_Project_bc_mut_gatk_hc_dbsnp2.tsv", header=TRUE)
sel = c("Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2", "dbSNP_RS", "Tumor_Sample_Barcode", "HGVSc", "HGVSp", "SYMBOL", "Hugo_Symbol", "Feature", "SIFT", "PolyPhen", "FILTER")
maf3 = maf2[maf2$SYMBOL %in% gatk$symbol, sel]

# filter common variants, af>0.01, gnomad
maf3 = maf3[maf3$FILTER!="common_variant", colnames(maf3)!="FILTER"]
# take out mutations occuring in more than 20% of samples
maf3$chr_start = paste(maf3$Chromosome, maf3$Start_Position) 
maf3 = maf3[order(maf3$chr_start),]
tmp = as.data.frame(table(maf3$chr_start))
tmp = as.vector(tmp$Var1[tmp$Freq>ceiling(nrow(samples)/5)])
maf3 = maf3[!maf3$chr_start %in% tmp, colnames(maf3)!="chr_start"]

pdf(paste('plots/waterfall',expname,'.pdf',sep=''), height=14,width=12)
waterfall(maf3, maxGenes = 50, mainDropMut=T, plotMutBurden = FALSE, mainXlabel=TRUE, section_heights=c(0.1,10))
dev.off()

write.table(maf3, paste('tables/waterfall',expname,'.csv',sep=''), row.names=F, quote=F, sep="\t")

## transform csv files to xls tables
setwd(paste(getwd(),'/tables',sep=''))
system('../csv2xlsNGS.pl')
system('rm *.csv')
setwd(gsub('/tables','',getwd()))


# distribution of mutation types
maf3 = readxl::read_excel("tables/waterfall_Project_bc_mut_gatk_hc_dbsnp_all_genes.xls")
pdf(paste("plots/barplot_mut_types_waterfall", expname, ".pdf", sep = ""), width = 3, height = 5)
tab = maf3 %>% count(Variant_Classification) %>% as.data.frame()
p1 <- ggplot(data = tab, aes(x = Variant_Classification, y = n)) +
    geom_bar(stat = "identity", position=position_dodge()) + #geom_jitter(width = 0.2) +
    theme_bw(base_size = 12, base_family = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Mutation types") +
    xlab("") + ylab("Number of mutations")
p1
dev.off()
pdf(paste("plots/pie_mut_types_waterfall", expname, ".pdf", sep = ""), width = 5, height = 3)
ggplot(tab, aes(x="", y=n, fill=Variant_Classification)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()  +
  scale_fill_brewer(palette="Paired")
dev.off()


pdf(paste("plots/barplot_mut_types_samples_waterfall", expname, ".pdf", sep = ""), width = 5, height = 10)
tab = maf3 %>% count(Variant_Classification) %>% as.data.frame()
p1 <- ggplot(data = tab, aes(x = Variant_Classification, y = n)) +
    geom_bar(stat = "identity", position=position_dodge()) + #geom_jitter(width = 0.2) +
    theme_bw(base_size = 12, base_family = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    ggtitle("Mutation types") +
    xlab("") + ylab("Number of mutations")
tab = maf3 %>% count(Tumor_Sample_Barcode) %>% as.data.frame()
p2 <- ggplot(data = tab, aes(x = Tumor_Sample_Barcode, y = n)) +
    geom_bar(stat = "identity", position=position_dodge()) + #geom_jitter(width = 0.2) +
    # scale_y_log10(limits=c(10^4,10^6)) +
    theme_bw(base_size = 12, base_family = "") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    # facet_wrap(~filter, nrow = 3) + 
    # scale_fill_brewer(palette = "Paired") +
    ggtitle("") +
    xlab("") + ylab("Number of mutations")
grid.arrange(p2, p1)
dev.off()


