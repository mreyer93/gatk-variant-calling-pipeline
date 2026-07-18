# differential expression for the MPL dataset
library("DESeq2")
library("BiocParallel")
library("apeglm")
library("ggplot2")
library("vsn")
library("reshape2")
library("ashr")
library("biomaRt")
library(sva)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
# setup parallelization
register(MulticoreParam(4))
parallel <- TRUE

count.f <- "~/pcloud_sync/project/rnaseq/mpl_experiment/results/counts/all.tsv"
sample.f <- "~/pcloud_sync/project/rnaseq/mpl_experiment/config/samples.tsv"
units.f <- "~/pcloud_sync/project/rnaseq/mpl_experiment/config/units.tsv"
outdir <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex'

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(count.f, header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(sample.f, header=TRUE, row.names="sample_name", check.names=FALSE)
nsamp <- nrow(coldata)
coldata <- coldata[colnames(cts)[1:nsamp], ,drop=F]
all(rownames(coldata) == colnames(cts)[1:nsamp])

# change rownames of count table first
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene.ids <- rownames(cts)
foo <- getBM(attributes=c('ensembl_gene_id',
                          'external_transcript_name'),
             filters = 'ensembl_transcript_id',
             values = gene.ids,
             mart = ensembl)
foo.gene <- getBM(attributes=c('ensembl_gene_id',
                               'external_gene_name'),
                  filters = 'ensembl_gene_id',
                  values = gene.ids,
                  mart = ensembl)
foo.desc <- getBM(attributes=c('ensembl_gene_id',
                               'description'),
                  filters = 'ensembl_gene_id',
                  values = gene.ids,
                  mart = ensembl)

ens.to.gene <- setNames(foo.gene$external_gene_name, foo.gene$ensembl_gene_id)
ens.to.desc <- setNames(foo.desc$description, foo.gene$ensembl_gene_id)
gene.to.desc <- setNames(ens.to.desc, ens.to.gene[names(ens.to.desc)])

cts$gene <- ens.to.gene[rownames(cts)]
dup.genes <- unique(cts$gene[duplicated(cts$gene)])
# make empty strings in to NA with a number
nna <- sum(cts$gene=='')
cts$gene[cts$gene==''] <- paste('NA', 1:nna, sep='_')
# aggregate counts across the same gene
cts.agg <- aggregate.data.frame(cts[,1:nsamp], by=list(cts$gene), FUN= sum)
rownames(cts.agg) <- unique(cts.agg$Group.1)
cts.agg <- cts.agg[,2:(nsamp+1)]

dds <- DESeqDataSetFromMatrix(countData=cts.agg,
                              colData=coldata,
                              design= ~ condition)
# Set reference factor levels
dds$condition <- relevel(dds$condition, ref = "untrt")
# dds$genotype <- relevel(dds$genotype, ref = "wt")

# only keep rows >= 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# DO IT
dds <- DESeq(dds, parallel=parallel)
resultsNames(dds)

# PCA PLOT # WITH ALL SAMPLES
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "genotype", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= batch, shape=interaction(genotype, condition))) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='MPL experiment plot: non batch corrected data')

# PCA PLOT # WITH BATCH1 SAMPLES
dds <- DESeqDataSetFromMatrix(countData=cts.agg[,coldata$batch=='SRA'],
                              colData=coldata[coldata$batch=='SRA',],
                              design= ~ genotype)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "genotype", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData,  aes(PC1, PC2, color= batch, shape=interaction(genotype, condition))) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='MPL experiment plot: non batch corrected data')

# PCA PLOT # WITH BATCH2 SAMPLES
dds <- DESeqDataSetFromMatrix(countData=cts.agg[,coldata$batch=='BOX'],
                              colData=coldata[coldata$batch=='BOX',],
                              design= ~ condition)
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "genotype", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= batch, shape=interaction(genotype, condition))) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='MPL experiment plot: non batch corrected data')



# doing COMBAT for normalization
# adjusted_counts <- ComBat_seq(as.matrix(cts.agg), batch=coldata$batch, group=NULL, covar_mod=coldata[,c('condition', 'genotype')])
adjusted_counts <- ComBat_seq(as.matrix(cts.agg), batch=coldata$batch)
# and then re-do deseq
dds <- DESeqDataSetFromMatrix(countData=adjusted_counts,
                              colData=coldata,
                              design= ~ condition + genotype)
# Set reference factor levels
dds$condition <- relevel(dds$condition, ref = "untrt")
dds$genotype <- relevel(dds$genotype, ref = "wt")

# only keep rows >= 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# DO IT
dds <- DESeq(dds, parallel=parallel)
resultsNames(dds)

# PCA PLOT
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "genotype", "batch"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= batch, shape=interaction(genotype, condition))) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='MPL experiment plot: BATCH CORRECTED')

# raw pca on scaled data
a <- prcomp(scale(t(counts(dds, normalized=T))))
ap <- cbind(data.frame(PC1=a$x[,1], PC2=a$x[,2]), coldata)
ggplot(ap, aes(PC1, PC2, color= batch, shape=interaction(genotype, condition))) +
    geom_point(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='MPL experiment plot: BATCH CORRECTED')


# dispersion plot
plotDispEsts(dds)


test.condition.treat <- 'condition_trt_vs_untrt'
test.condition.genotype <- 'genotype_MPLW515L_vs_wt'
# test.condition.interaction <- 'conditionLSDi.genotypeJAK2V617F'


###############################################################################
# TREATMENT EFFECT ON BATCH CORRECTED #########################################
###############################################################################
test.contrast <- c(0,1,0)
res.treat <- results(dds, contrast = test.contrast, alpha = 0.01,)
resOrdered.treat <- res.treat[order(res.treat$pvalue),]
resLFC.treat <- lfcShrink(dds, contrast =  test.contrast, type="ashr", res = res.treat)
resLFCOrdered.treat <- resLFC.treat[order(resLFC.treat$pvalue),]
summary(resLFC.treat)
plotMA(resLFC.treat, ylim=c(-3,3), main="LSDi effect")

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 9
plot.genes <- counts(dds)[order(resLFC.treat$padj)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$genotype <- coldata[pgm$sample,"genotype"]
pgm$condition <- factor(coldata[pgm$sample,"condition"], levels=c('untrt', 'trt'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))
pgm$batch <- factor(coldata[pgm$sample,"batch"])

ggplot(pgm, aes(x=interaction(condition,genotype, batch), y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))


# this doesn't look good  - much of the signal is still in batch.
# will move to just quantifying effect within each batch so that we can
# get some reliable results
outf.treat <- file.path(outdir, 'diffex_treatment.tsv')
out.treat <- data.frame(name = rownames(resLFCOrdered.treat),
                        description = gene.to.desc[rownames(resLFCOrdered.treat)],
                        signif(as.data.frame(resLFCOrdered.treat),5))

write.table(out.treat, outf.treat, sep='\t', quote=F, col.names = T, row.names = F)

###############################################################################
# UNCORRECTED FOR EACH BATCH SEPARATE #########################################
###############################################################################
# PCA PLOT # WITH BATCH1 SAMPLES
dds.sra <- DESeqDataSetFromMatrix(countData=cts.agg[,coldata$batch=='SRA'],
                              colData=coldata[coldata$batch=='SRA',],
                              design= ~ genotype)

# PCA PLOT # WITH BATCH2 SAMPLES
dds.box <- DESeqDataSetFromMatrix(countData=cts.agg[,coldata$batch=='BOX'],
                              colData=coldata[coldata$batch=='BOX',],
                              design= ~ condition)

# Set reference factor levels
dds.box$condition <- relevel(dds.box$condition, ref = "untrt")
dds.sra$genotype <- relevel(dds.sra$genotype, ref = "wt")

# only keep rows >= 10 counts
keep.sra <- rowSums(counts(dds.sra)) >= 10
keep.box <- rowSums(counts(dds.box)) >= 10
dds.sra <- dds.sra[keep.sra,]
dds.box <- dds.box[keep.box,]

# DO IT
dds.sra <- DESeq(dds.sra, parallel=parallel)
dds.box <- DESeq(dds.box, parallel=parallel)
resultsNames(dds.sra)
resultsNames(dds.box)

# export counts table
counts.export.sra <- data.frame(name=rownames(counts(dds.sra)),
                            description =gene.to.desc[rownames(counts(dds.sra))],
                            round(counts(dds.sra, normalized=T), 3))
outf.counts.sra <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/counts/normalized_counts_export_SRA.tsv'
write.table(counts.export.sra, outf.counts.sra, sep='\t', quote=F, row.names = F, col.names = T)
counts.export.box <- data.frame(name=rownames(counts(dds.box)),
                            description =gene.to.desc[rownames(counts(dds.box))],
                            round(counts(dds.box, normalized=T), 3))
outf.counts.box <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/counts/normalized_counts_export_BOX.tsv'
write.table(counts.export.box, outf.counts.box, sep='\t', quote=F, row.names = F, col.names = T)

# calculate means of each condition
coldata.sra <- coldata[coldata$batch == 'SRA',]
coldata.sra$final.group <- paste(coldata.sra$genotype, coldata.sra$condition, sep='_')
norm.counts.sra <- counts(dds.sra, normalized=T)
mean.group.counts.sra <- do.call(cbind, lapply(unique(coldata.sra$final.group), function(x) {
    rowMeans(norm.counts.sra[,rownames(coldata.sra)[which(coldata.sra$final.group==x)]])}))
colnames(mean.group.counts.sra) <- unique(coldata.sra$final.group)
counts.mean.export.sra <- data.frame(name=rownames(mean.group.counts.sra),
                                 description =gene.to.desc[rownames(mean.group.counts.sra)],
                                 dispersion=round(dispersions(dds.sra),4),
                                 round(mean.group.counts.sra, 3))
outf.counts.mean.sra <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/counts/normalized_counts_group_mean_SRA.tsv'
write.table(counts.mean.export.sra, outf.counts.mean.sra, sep='\t', quote=F, row.names = F, col.names = T)
# for BOX
coldata.box <- coldata[coldata$batch == 'BOX',]
coldata.box$final.group <- paste(coldata.box$genotype, coldata.box$condition, sep='_')
norm.counts.box <- counts(dds.box, normalized=T)
mean.group.counts.box <- do.call(cbind, lapply(unique(coldata.box$final.group), function(x) {
    rowMeans(norm.counts.box[,rownames(coldata.box)[which(coldata.box$final.group==x)]])}))
colnames(mean.group.counts.box) <- unique(coldata.box$final.group)
counts.mean.export.box <- data.frame(name=rownames(mean.group.counts.box),
                                 description =gene.to.desc[rownames(mean.group.counts.box)],
                                 dispersion=round(dispersions(dds.box),4),
                                 round(mean.group.counts.box, 3))
outf.counts.mean.box <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/counts/normalized_counts_group_mean_BOX.tsv'
write.table(counts.mean.export.box, outf.counts.mean.box, sep='\t', quote=F, row.names = F, col.names = T)

###############################################################################
# Test genotype effect - SRA samples only
contrast.genotype <- c(0,1)
res.genotype <- results(dds.sra, contrast = contrast.genotype, alpha = 0.01,)
resOrdered.genotype <- res.genotype[order(res.genotype$pvalue),]
resLFC.genotype <- lfcShrink(dds.sra, coef= test.condition.genotype, type="ashr", res = res.genotype)
resLFCOrdered.genotype <- resLFC.genotype[order(resLFC.genotype$pvalue),]
summary(resLFC.genotype)
plotMA(resLFC.genotype, ylim=c(-3,3), main="Genotype effect: SRA samples")

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 9
plot.genes <- counts(dds.sra)[order(resLFC.genotype$padj)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$genotype <- coldata[pgm$sample,"genotype"]
pgm$condition <- factor(coldata[pgm$sample,"condition"], levels=c('wt', 'MPLW515L'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))
ggplot(pgm, aes(x=genotype, y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75)) +
    labs(main='Differential genes: SRA samples')

outf.genotype <- file.path(outdir, 'diffex_genotype_SRA.tsv')
out.genotype <- data.frame(name = rownames(resLFCOrdered.genotype),
                           description = gene.to.desc[rownames(resLFCOrdered.genotype)],
                           signif(as.data.frame(resLFCOrdered.genotype),5))

write.table(out.genotype, outf.genotype, sep='\t', quote=F, col.names = T, row.names = F)

###############################################################################
# Test treatment effect - BOX mutant samples only
contrast.treat <- c(0,1)
res.treat <- results(dds.box, contrast = contrast.treat, alpha = 0.01)
resOrdered.treat <- res.treat[order(res.treat$pvalue),]
resLFC.treat <- lfcShrink(dds.box, contrast =  test.contrast, type="ashr", res = res.treat)
resLFCOrdered.treat <- resLFC.treat[order(resLFC.treat$pvalue),]
summary(resLFC.treat)
plotMA(resLFC.treat, ylim=c(-3,3), main="Treatment effect: BOX samples (MPLW515L)")

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 9
plot.genes <- counts(dds.box)[order(resLFC.treat$padj)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$genotype <- coldata[pgm$sample,"genotype"]
pgm$condition <- factor(coldata[pgm$sample,"condition"], levels=c('untrt', 'trt'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))
pgm$batch <- factor(coldata[pgm$sample,"batch"])

ggplot(pgm, aes(x=interaction(genotype, condition), y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75)) +
    labs(main='Differential genes: BOX samples')



outf.treat <- file.path(outdir, 'diffex_treatment_BOX.tsv')
out.treat <- data.frame(name = rownames(resLFCOrdered.treat),
                        description = gene.to.desc[rownames(resLFCOrdered.treat)],
                        signif(as.data.frame(resLFCOrdered.treat),5))

write.table(out.treat, outf.treat, sep='\t', quote=F, col.names = T, row.names = F)
