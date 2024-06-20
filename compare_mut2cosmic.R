###################################################
## compare mutations with cosmic
###################################################

library(psych)
library(tidyr)
library(ggplot2)
library(data.table)

expname = paste("_",basename(getwd()), "_gatk_hc_dbsnp2", sep="")


# cell line samples
samples = read.csv("file2sample.csv")
samples = samples[!duplicated(samples$sampleName), ]


# A tab separated table of all the point mutations in the Cosmic Cell Lines Project from the current release.
tab = read.csv(gzfile("../data/CosmicCLP_MutantExport_v97_2023_01.tsv.gz"), sep="\t")
colnames(tab) = gsub("\\.", "_", colnames(tab))
tab$Sample_name = toupper(tab$Sample_name)
sum(samples$sampleName %in% tab$Sample_name)
sort(samples$sampleName[(!samples$sampleName %in% tab$Sample_name)])
# change cell lines names to dsmz notation
tab$Sample_name[tab$Sample_name == "EoL-1-cell" ] = "EOL-1"
tab$Sample_name[tab$Sample_name == "KU812" ] = "KU-812"
tab$Sample_name[tab$Sample_name == "OCI-LY-19" ] = "OCI-LY19"
tab$Sample_name[tab$Sample_name == "Set2" ] = "SET-2"
tab$Sample_name[tab$Sample_name == "HCC1599" ] = "HCC-1599"
tab$Sample_name[tab$Sample_name == "HCC1937" ] = "HCC-1937"
tab$Sample_name[tab$Sample_name == "HS-578-T" ] = "HS-578T"
tab$Sample_name[tab$Sample_name == "MCF7" ] = "MCF-7"
tab$Sample_name[tab$Sample_name == "T47D" ] = "T-47D"
sum(samples$sampleName %in% tab$Sample_name)
sort(samples$sampleName[(!samples$sampleName %in% tab$Sample_name)])
# focus on specific cancer cell lines
tab = tab[tab$Sample_name %in% samples$sampleName,]
# get chromosomal location
tmp = t(sapply(tab$Mutation_genome_position, function(x) strsplit(x,":")[[1]]))
tmp[,1] = paste("chr", tmp[,1], sep="")
tmp[,2] = sapply(tmp[,2], function(x) strsplit(x, "-")[[1]][[1]])
tab$Chrom = tmp[,1]
tab$Pos = as.numeric(tmp[,2])

# focus on cosmic verified mutations + DSMZ (Braunschweig, Germany) + without unknown mutations => 10 cell lines remaining
tab = tab[tab$Mutation_verification_status == "Verified",]
tab = tab[tab$Institute_Address == " Braunschweig, Germany",]
tab = tab[tab$Mutation_Description != "Unknown",]
tab = tab[!duplicated(paste(tab$Sample_name, tab$Chrom, tab$Pos, tab$HGVSG)), ]

mut = tab[,c(1:2, 5, 13, 21:23, 26:27, 36:42)]
mut$symbol = sapply(mut$Gene_name, function(x) strsplit(x, "_")[[1]][1])

# get vcf and annotation files
files = paste("snp/", sort(unique(mut$Sample_name)), sep="")
files_vep = paste(files, ".hc.pass2.lcr.dbsnp.vcf.maf.gz", sep="") # maf format
files = paste(files, ".hc.pass2.lcr.dbsnp.vcf.gz", sep="")
name = sort(unique(mut$Sample_name))
minCount = 5 # minimal depth
# add info, whether mutation is found by gatk, first set "no"
mut$gatk= "no"

# all variants of the cell lines in one table
for (i in seq_along(name)) {
    print(paste(files[i], date()))
    ## gatk results
    gk = read.csv(files[i], header = T, stringsAsFactors = FALSE, skip = 228, sep = "\t")
    # all strand bias, SNPcluster etc filtered by gatk
    gk = gk[gk$FILTER=="PASS", ]
    names(gk)[1] = "CHROM"
    gk$AD = sapply(gk[, 10], function(x) strsplit(x, ":")[[1]][2])
    gk$DP = sapply(gk$AD, function(x) sum(as.numeric(strsplit(x, ",")[[1]])))
    gk = gk[gk$DP >= minCount, ]
    gk$AD = sapply(gk$AD, function(x) as.numeric(strsplit(x, ",")[[1]][2]))
    colnames(gk)[(ncol(gk)-1):ncol(gk)] = paste(name[i],c("AD","DP"), sep="_")
    gk = gk[, -c(3, 6:10)]
    # create id for overlapping with vep; caution with deletion/insertion
	gk$id = paste(gk$CHROM, gk$POS, sep = "_")

    ## load variant effect information (vep), maf format
    # take out mutations gnomad af>0.01 => in maf FILTER="common_variant"
    # + filter coding variants
    vep = read.csv(gzfile(files_vep[i]), header=T, stringsAsFactors=F, sep = "\t", comment.char="#")
    vep = vep[grep("protein_coding", vep$BIOTYPE), c(1:16, 35:43, 47:75, 77:87, 93:114)]
    vep = vep[vep$HGVSp!="", ]
    vep = vep[vep$FILTER != "common_variant", ] # gnomad filter
    # => --af => phase 3 1000 genomes
    vep = vep[is.na(vep$AF) | vep$AF<0.01, ]
    # keep variants only with minimum depth via gatk data
    vep$id = paste(vep$Chromosome, vep$vcf_pos, sep="_") # vep$vcf_pos: already -1 for deletions
    vep = vep[vep$id %in% gk$id,]

    # overlap to cosmic
    mut$gatk[mut$Chrom %in% vep$Chromosome & mut$Mutation_AA %in% vep$HGVSp_Short & mut$Sample_name == name[i]] = "yes"
    tmp = mut[mut$Sample_name == name[i], ]
}

mut = mut[order(mut$Sample_name), ]

# cancer gene census, curated cancer genes
cgc = read.csv(gzfile("../data/Cosmic_MutantCensus_v98_GRCh38.tsv.gz"), sep="\t")
cgc = cgc[cgc$SAMPLE_NAME %in% unique(mut$Sample_name),] # COLO-824 not found
sort(unique(cgc$MUTATION_AA))
tmp = match(mut$Mutation_AA, cgc$MUTATION_AA, nomatch=0)
tmp = tmp[tmp!=0]

mut$cgc = "no"
for (i in 1:nrow(mut)) {
    tmp = cgc[cgc$MUTATION_AA %in% mut$Mutation_AA[i] & cgc$SAMPLE_NAME %in% mut$Sample_name & cgc$GENE_SYMBOL %in% mut$symbol,]
    if (nrow(tmp>0)) {
        mut$cgc[i] = "yes"
    }
}

sum(mut$gatk == "no" & mut$cgc == "yes")
sum(mut$gatk == "yes" & mut$cgc == "yes")
sum(mut$cgc == "yes")
#  => 10 in mut table, all 10 found by gatk

mut2 = mut[!duplicated(paste(mut$Sample_name, mut$HGVSG)), c(3,5,6, 13:19)]


# check, if missed mutants are expressed
# tpm better than nexprs => all expr values within one sample add to 1M counts; this different to nexprs
tpm = read.csv("../tables/TPM_salmon_Project_bc.tsv", sep = "\t")
colnames(tpm) = gsub("\\.", "-", colnames(tpm))
mut2$symbol[!mut2$symbol %in% tpm$external_gene_name]
# "FAM129B" "TMEM206" => current names: PACC1, NIBAN2
mut2$symbol[mut2$symbol == "FAM129B"] = "PACC1"
mut2$symbol[mut2$symbol == "TMEM206"] = "NIBAN2"
mut2 = mut2[order(mut2$Pos),]
mut2 = mut2[order(mut2$Chrom),]
mut2 = mut2[order(mut2$Sample_name),]
genes = unique(mut2$symbol)


# non-expressed genes in variant table
mut2$expr = "<=1 tpm"
for (i in 1:nrow(mut2)) {
    tmp = tpm[tpm$external_gene_name %in% mut2$symbol[i], colnames(tpm) %in% mut2$Sample_name[i]]
    if (tmp>1) {
        mut2$expr[i] = ">1 tpm"
    }
}

# Missing variants in filtered rna edit sites?
redit = read.csv("../data/TABLE1_hg38_chr_pos.txt", sep="\t", header=FALSE)
mut2$redit = "no"
for (i in 1:nrow(mut2)) {
    tmp = redit[redit$V1 %in% mut2$Chrom[i] & redit$V2 %in% mut2$Pos[i],]
    if(i%%50 ==0) print(paste(date(), i))
    if (nrow(tmp)>0) {
        mut2$redit[i] = "yes"
    }
}
mut2[mut2$redit=="yes",]
rm(redit)

# Missing variants in filtered low complexity regions?
lcr = read.csv(gzfile("../data/LCR-hs38_with_chr.bed.gz"), sep="\t", header=FALSE)
mut2$LCR = "out"
for (i in 1:nrow(mut2)) {
    tmp = lcr[lcr$V1 %in% mut2$Chrom[i] & lcr$V2 <= mut2$Pos[i] & lcr$V3 > mut2$Pos[i],]
    if(i%%50 ==0) print(paste(date(), i))
    if (nrow(tmp)>0) {
        mut2$LCR[i] = "in"
    }
}
mut2[mut2$LCR=="in",]
rm(lcr)

# Missing variants in filtered common dbsnp sites?
dbsnp = read.csv("../data/GCF_000001405.40.common.chr.tsv", sep="\t", header=FALSE)
mut2$dbsnp = "no"
for (i in 1:nrow(mut2)) {
    tmp = dbsnp[dbsnp$V1 %in% mut2$Chrom[i] & dbsnp$V2 %in% mut2$Pos[i],]
    if(i%%50 ==0) print(paste(date(), i))
    if (nrow(tmp)>0) {
        mut2$dbsnp[i] = "yes"
    }
}
mut2[mut2$dbsnp=="yes",]
rm(dbsnp)

# Missing variants in filtered common 1kG sites?
kG = read.csv("../data/1000G_phase3_v4_20130502.sites.hg38.af01.vcf", sep="\t", skip=3527)
mut2$kG = "no"
for (i in 1:nrow(mut2)) {
    tmp = kG[kG$V1 %in% mut2$Chrom[i] & kG$V2 %in% mut2$Pos[i],]
    if(i%%50 ==0) print(paste(date(), i))
    if (nrow(tmp)>0) {
        mut2$kG[i] = "yes"
    }
}
mut2[mut2$kG=="yes",]
rm(kG)

# Missing variants in filtered >20% samples sites?
gatk = readxl::read_excel("tables/snv_indel_maf_less20_Project_bc_mut_gatk_hc_dbsnp.tsv")
mut2$Freq = ""
for (i in 1:nrow(mut2)) {
    if(mut2$gatk[i]=="no") next
    tmp = gatk[gatk$Chromosome %in% mut2$Chrom[i] & 
        gatk$sampleID %in% mut2$Sample_name[i] &
        (gatk$Start_Position %in% mut2$Pos[i] | gatk$HGVSp_Short == mut2$Mutation_AA[i]),]
    if(i%%50 ==0) print(paste(date(), i))
    if (nrow(tmp)>0) {
        mut2$Freq[i] = "<20%"
    }
}
rm(gatk)


# Depth at mutation sites via samtools, preparation of needed bed file
# as bed file for samtools for getting depth at each mutation sites
tmp = mut2[,c("Chrom", "Pos", "Pos")]
write.table(tmp, paste("tables/cosmic_mut", expname, ".bed", sep = ""),
    row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE)
# as bed file for samtools for getting depth at each mutation sites, 0 based
tmp$Pos = tmp$Pos - 1
write.table(tmp, paste("tables/cosmic_mut", expname, ".0based.bed", sep = ""),
    row.names = FALSE, col.names = FALSE, sep = "\t", quote=FALSE)
# perform samtools bedcov on server and proceed further afterwards
write.table(mut2, paste("tables/cosmic_mut", expname, ".tsv", sep = ""),
    row.names = FALSE, sep = "\t")

# save cosmic mutation in bc mutation table
tab = read.csv(paste('tables/snv_indel_maf_Project_bc_mut_gatk_hc_dbsnp.tsv',sep=''), header=TRUE, sep="\t")
tab$COSMIC = FALSE
for (i in 1:nrow(mut2)) {
    tmp = tab[tab$sampleID == mut2$Sample_name[i] & tab$Chromosome == mut2$Chrom[i] & tab$HGVSp_Short == mut2$Mutation_AA[i],]
    if(nrow(tmp)>0) tab$COSMIC[tab$sampleID == mut2$Sample_name[i] & tab$Chromosome == mut2$Chrom[i] & tab$HGVSp_Short == mut2$Mutation_AA[i]] = TRUE
}
write.table(tab, paste("tables/snv_indel_maf_Project_bc_mut_gatk_hc_dbsnp.tsv", sep = ""),
    row.names = FALSE, sep = "\t")


sessionInfo()

