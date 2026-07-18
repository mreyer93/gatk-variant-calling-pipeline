#  clustering of patient VCFs to get a sense of tumor-normal relationships
# whole exome and targeted EGA data
library(vcfR)
library(proxy)
library(SNPRelate)
library(SeqArray)
library(ggplot2)
library(ggdendro)

vcf_file <- '~/data/EGA_data/02_variants_reference/merged_filtered.vcf.gz'
vcf_file_no0 <- '~/data/EGA_data/02_variants_reference/merged_filtered_no0.vcf.gz'
vcf <- read.vcfR( vcf_file, verbose = FALSE )
vcf.no0 <- read.vcfR( vcf_file_no0, verbose = FALSE )


gt <- extract.gt(vcf)
dp <- extract.gt(vcf,element = "DP" ,as.numeric = T)
dp.rs <- rowSums(!is.na(dp), na.rm = T)
dp.cs <- colSums(!is.na(dp), na.rm = T)
hist(dp.rs)
plot(density(dp.rs))
# only keep variants where we have >2 samples?
keep.variants <- which(dp.rs >2)
vcf <- vcf[keep.variants, ]


gt <- extract.gt(vcf)
gt.table <- table(gt)
# change the phased / unphased
gt <- gsub('/', '|', gt)
gt.table <- table(gt)
sort(gt.table, decreasing=T)

gt.jd <- dist(gt, method='manhattan', by_rows = F, pairwise = F)
gt.jd.mat <- as.matrix(gt.jd)
gt.jd.mat[1:5,1:5]

# doing with snprelate instead
vcf_file <- '~/data/EGA_data/02_variants_reference/merged_filtered.vcf.gz'
vcf_file_no0 <- '~/data/EGA_data/02_variants_reference/merged_filtered_no0.vcf.gz'
seqVCF2GDS(vcf_file, "test.gds",parallel = 64)

genofile <- seqOpen("test.gds")
gt <- seqGetData(genofile, "genotype")

# exome data on cluster
vcf_file <- '~/data/EGA_data/02_variants_reference_exome/merged_filtered.vcf.gz'
vcf_file_no0 <- '~/data/EGA_data/02_variants_reference_exome/merged_filtered_no0.vcf.gz'
seqVCF2GDS(vcf_file, "test.gds", parallel = 64)
seqVCF2GDS(vcf_file_no0, "test_no0.gds",parallel = 64)
genofile <- seqOpen("test.gds")
genofile_no0 <- seqOpen("test_no0.gds")
# gt <- seqGetData(genofile, "genotype")


snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2, num.thread= 64)
snpset.id <- unlist(unname(snpset))
saveRDS(snpset.id, 'snpsed.id.rds')
length(snpset.id)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=32)
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

pca <- ggplot(tab, aes(x=EV1, y=EV2, label=sample.id)) +
    geom_point() +
    geom_text(nudge_y = 0.05)

remove.samples <- c('PD6637a2', 'PD6647a2')
tab.removed <- tab[!(tab$sample.id %in% remove.samples), ]

ggplot(tab.removed, aes(x=EV1, y=EV2, label=sample.id)) +
    geom_point() +
    geom_text(nudge_y = 0.005)

# distance with thier functions
ds <- snpgdsDiss(genofile, snp.id=snpset.id)
ds$diss[1:5,1:5]

ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, snp.id = snpset.id, num.thread=64))
rv <- snpgdsCutTree(ibs.hc)

pdf('~/pcloud_sync/project/sample_vcf_clustering_trees/EGA_IBS_tree.pdf', height=45, width=6)
plot.new()
p <- ggdendrogram(rv$dendrogram, drop=T, rotate = T)
print(p)
dev.off()

snpset_no0 <- snpgdsLDpruning(genofile_no0, ld.threshold=0.2, num.thread= 64)
snpset.id_no0 <- unlist(unname(snpset_no0))
saveRDS(snpset.id_no0, 'snpsed.id_no0.rds')
pca_no0 <- snpgdsPCA(genofile_no0, snp.id=snpset.id_no0, num.thread=32)
tab_no0 <- data.frame(sample.id = pca_no0$sample.id,
                  EV1 = pca_no0$eigenvect[,1],    # the first eigenvector
                  EV2 = pca_no0$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

pca_no0 <- ggplot(tab_no0, aes(x=EV1, y=EV2, label=sample.id)) +
    geom_point() +
    geom_text(nudge_y = 0.05)

pdf('test_pca_no0.pdf')
print(pca_no0)


# distance with thier functions
ds_no0 <- snpgdsDiss(genofile_no0, snp.id=snpset.id_no0)

ibs_no0 <- snpgdsIBS(genofile_no0, snp.id = snpset.id_no0, num.thread=64)
ibs.hc_no0 <- snpgdsHCluster(ibs_no0)
rv_no0 <- snpgdsCutTree(ibs.hc_no0)

pdf('test_tree_no0.pdf', height=45, width=6)
plot.new()
p <- ggdendrogram(rv_no0$dendrogram, drop=T, rotate = T)
print(p)
dev.off()
