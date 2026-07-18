# processing and clustering CTP, etc data
library(SNPRelate)
library(SeqArray)
library(ggplot2)
library(ggdendro)

# CTP_102
seqClose(genofile)
vcf.fn <- '~/data/ctp_102_calling/output/02_variants_reference/filtered/merged/merged_filtered.vcf'
# CTP_201
seqClose(genofile)
vcf.fn <- '~/data/ctp_201_calling/output/02_variants_reference/filtered/merged/merged_filtered.vcf'
# EGA DATA
seqClose(genofile)
vcf.fn <- '~/data/EGA_data/02_variants_reference/filtered/merged/merged_filtered.vcf'

seqVCF2GDS(vcf.fn, "test2.gds")
genofile <- seqOpen("test2.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(unname(snpset))
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=4)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

# get raw genotype information
gg <- snpgdsGetGeno(genofile)
samples <- seqGetData(genofile, 'sample.id')
rownames(gg) <- samples
# plot(hclust(dist(gg)))
x
model <- hclust(dist(gg, method = 'manhattan'), "complete")
# pdf('~/ctp102_tree.pdf', height=35, width=6)
pdf('~/EGA_tree.pdf', height=45, width=6)
# pdf('~/ctp201_tree.pdf', height=8, width=6)
ggdendrogram(model, rotate = TRUE, size = 2)
dev.off()

# cut tree height
h=35
pt.data <- data.frame(sample=samples, pt=sapply(samples, function(x) strsplit(x, split='_')[[1]][1]), cluster=cutree(model, h=h))

# are there any clusters where different pts are in it?
for (clus in unique(pt.data$cluster)){
    pt.data.sub <- pt.data[pt.data$cluster==clus,]
    if(!all(pt.data.sub$pt == pt.data.sub$pt[1])){
        print(pt.data.sub)
    }
}

# 020-102 and 020-103 cluster very close together - what's up with that?
# 009-105_Dx_BS clusters with the rest of 021-101
