# processing count data into DESeq2
# count data generated with STAR alignment via the snakemake pipeline
# using custom analysis script for more flexibility
library("DESeq2")
library("BiocParallel")
library("apeglm")
library("ggplot2")
library("vsn")
library("reshape2")
library("ashr")
library("biomaRt")
httr::set_config(httr::config(ssl_verifypeer = FALSE))
# setup parallelization
register(MulticoreParam(4))
parallel <- TRUE

count.f <- "~/pcloud_sync/project/rnaseq/batch1/results/counts/all.tsv"
sample.f <- "~/pcloud_sync/project/rnaseq/batch1/config/samples.tsv"
units.f <- "~/pcloud_sync/project/rnaseq/batch1/config/units.tsv"
outdir <- '~/pcloud_sync/project/rnaseq/batch1/results/diffex/ensembl_genes/'

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(count.f, header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(sample.f, header=TRUE, row.names="sample_name", check.names=FALSE)
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design= ~ condition * genotype)
# Set reference factor levels
dds$condition <- relevel(dds$condition, ref = "Vehicle")
dds$genotype <- relevel(dds$genotype, ref = "Bl6_wt")

# only keep rows >= 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# DO IT
dds <- DESeq(dds, parallel=parallel)
resultsNames(dds)

# PCA PLOT
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=genotype)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='Mouse RNA-seq PCA plot')

# dispersion plot
plotDispEsts(dds)

test.condition.treat <- 'condition_LSDi_vs_Vehicle'
test.condition.genotype <- 'genotype_JAK2V617F_vs_Bl6_wt'
test.condition.interaction <- 'conditionLSDi.genotypeJAK2V617F'

# export counts table
counts.export <- data.frame(name=rownames(counts(dds)),
                            round(counts(dds, normalized=T), 3))
outf.counts <- '~/pcloud_sync/project/rnaseq/batch1/results/counts/normalized_counts_export_ensembl.tsv'
write.table(counts.export, outf.counts, sep='\t', quote=F, row.names = F, col.names = T)

# calculate means of each condition
coldata$final.group <- paste(coldata$genotype, coldata$condition, sep='_')

norm.counts <- counts(dds, normalized=T)
mean.group.counts <- do.call(cbind, lapply(unique(coldata$final.group), function(x) {
    rowMeans(norm.counts[,rownames(coldata)[which(coldata$final.group==x)]])}))
colnames(mean.group.counts) <- unique(coldata$final.group)
counts.mean.export <- data.frame(name=rownames(mean.group.counts),
                                 dispersion=round(dispersions(dds),4),
                                 round(mean.group.counts, 3))
outf.counts.mean <- '~/pcloud_sync/project/rnaseq/batch1/results/counts/normalized_counts_group_mean_ensembl.tsv'
write.table(counts.mean.export, outf.counts.mean, sep='\t', quote=F, row.names = F, col.names = T)


###############################################################################
# TREATMENT EFFECT ############################################################
###############################################################################
test.contrast <- c(0,1,0,0.5)
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
pgm$condition <- factor(coldata[pgm$sample,"condition"], levels=c('Vehicle', 'LSDi'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))

ggplot(pgm, aes(x=interaction(condition,genotype), y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))

outf.treat <- file.path(outdir, 'diffex_treatment.tsv')
out.treat <- data.frame(name = rownames(resLFCOrdered.treat),
                        signif(as.data.frame(resLFCOrdered.treat),5))

write.table(out.treat, outf.treat, sep='\t', quote=F, col.names = T, row.names = F)

###############################################################################
# GENOTYPE EFFECT ############################################################
###############################################################################
test.contrast <- c(0,0,1,0.5)
res.genotype <- results(dds, contrast = test.contrast, alpha = 0.01,)
resOrdered.genotype <- res.genotype[order(res.genotype$pvalue),]
resLFC.genotype <- lfcShrink(dds, coef= test.condition.genotype, type="ashr", res = res.genotype)
resLFCOrdered.genotype <- resLFC.genotype[order(resLFC.genotype$pvalue),]
summary(resLFC.genotype)
plotMA(resLFC.genotype, ylim=c(-3,3), main="Genotype effect")

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 9
plot.genes <- counts(dds)[order(resLFC.genotype$padj)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$genotype <- coldata[pgm$sample,"genotype"]
pgm$condition <- factor(coldata[pgm$sample,"condition"], levels=c('Vehicle', 'LSDi'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))

ggplot(pgm, aes(x=interaction(condition,genotype), y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))

outf.genotype <- file.path(outdir, 'diffex_genotype.tsv')
out.genotype <- data.frame(name = rownames(resLFCOrdered.genotype),
                           signif(as.data.frame(resLFCOrdered.genotype),5))

write.table(out.genotype, outf.genotype, sep='\t', quote=F, col.names = T, row.names = F)

###############################################################################
# INTERACTION EFFECT ############################################################
###############################################################################
test.contrast <- c(0,0,0,1)
res.interaction <- results(dds, contrast = test.contrast, alpha = 0.01,)
resOrdered.interaction <- res.interaction[order(res.interaction$pvalue),]
resLFC.interaction <- lfcShrink(dds, contrast = test.contrast, type="ashr", res = res.interaction)
resLFCOrdered.interaction <- resLFC.interaction[order(resLFC.interaction$pvalue),]
summary(resLFC.interaction)
plotMA(resLFC.interaction, ylim=c(-3,3), main=test.condition.interaction)

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 20
plot.genes <- counts(dds)[rownames(resLFCOrdered.interaction)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$genotype <- coldata[pgm$sample,"genotype"]
pgm$condition <- factor(coldata[pgm$sample,"condition"], levels=c('Vehicle', 'LSDi'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))

ggplot(pgm, aes(x=interaction(condition,genotype), y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75))

outf.interaction <- file.path(outdir, 'diffex_interaction.tsv')
out.interaction <- data.frame(name = rownames(resLFCOrdered.interaction),
                              signif(as.data.frame(resLFCOrdered.interaction),5))

write.table(out.interaction, outf.interaction, sep='\t', quote=F, col.names = T, row.names = F)


###############################################################################
# NEW SPECIFIC COMPARISONS ####################################################
#  lsd1 in jak2 background ####################################################
#     lsd1 in wt mice      ####################################################
###############################################################################
# have to eliminate the other samples to do these exact comparisons
dds <- DESeqDataSetFromMatrix(countData=cts[, rownames(coldata)[coldata$genotype=='JAK2V617F']],
                              colData=coldata[coldata$genotype=='JAK2V617F',],
                              design= ~ condition)
# Set reference factor levels
dds$condition <- relevel(dds$condition, ref = "Vehicle")
# only keep rows >= 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# DO IT
dds <- DESeq(dds, parallel=parallel)
resultsNames(dds)
# PCA PLOT
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='Mouse RNA-seq PCA plot')

# dispersion plot
plotDispEsts(dds)

test.condition.treat <- 'condition_LSDi_vs_Vehicle'
res.treat <- results(dds, name=test.condition.treat, alpha = 0.01,)
resOrdered.treat <- res.treat[order(res.treat$pvalue),]
resLFC.treat <- lfcShrink(dds, type="ashr", res = res.treat)
resLFCOrdered.treat <- resLFC.treat[order(resLFC.treat$pvalue),]
summary(resLFC.treat)
plotMA(resLFC.treat, ylim=c(-3,3), main="LSDi effect in JAK2 Mutant")

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 9
plot.genes <- counts(dds)[order(resLFC.treat$padj)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$condition <- factor(coldata[pgm$sample,"condition"], levels=c('Vehicle', 'LSDi'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))

ggplot(pgm, aes(x=condition, y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75)) +
    labs(title = 'Treatment effect in JAK2 Mutant mice')

outf.treat <- file.path(outdir, 'diffex_treatment_JAK2V617F_ONY.tsv')
out.treat <- data.frame(name = rownames(resLFCOrdered.treat),
                        signif(as.data.frame(resLFCOrdered.treat),5))

write.table(out.treat, outf.treat, sep='\t', quote=F, col.names = T, row.names = F)

#     lsd1 in wt mice      ####################################################
###############################################################################
# have to eliminate the other samples to do these exact comparisons
dds <- DESeqDataSetFromMatrix(countData=cts[, rownames(coldata)[coldata$genotype=='Bl6_wt']],
                              colData=coldata[coldata$genotype=='Bl6_wt',],
                              design= ~ condition)
# Set reference factor levels
dds$condition <- relevel(dds$condition, ref = "Vehicle")
# only keep rows >= 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# DO IT
dds <- DESeq(dds, parallel=parallel)
resultsNames(dds)
# PCA PLOT
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='Mouse RNA-seq PCA plot')

# dispersion plot
plotDispEsts(dds)

test.condition.treat <- 'condition_LSDi_vs_Vehicle'
res.treat <- results(dds, name=test.condition.treat, alpha = 0.01,)
resOrdered.treat <- res.treat[order(res.treat$pvalue),]
resLFC.treat <- lfcShrink(dds, type="ashr", res = res.treat)
resLFCOrdered.treat <- resLFC.treat[order(resLFC.treat$pvalue),]
summary(resLFC.treat)
plotMA(resLFC.treat, ylim=c(-3,3), main="LSDi effect in Bl6_wt")

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 9
plot.genes <- counts(dds)[order(resLFC.treat$padj)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$condition <- factor(coldata[pgm$sample,"condition"], levels=c('Vehicle', 'LSDi'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))

ggplot(pgm, aes(x=condition, y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75)) +
    labs(title = 'Treatment effect in Bl6_wt mice')

outf.treat <- file.path(outdir, 'diffex_treatment_WT_ONY.tsv')
out.treat <- data.frame(name = rownames(resLFCOrdered.treat),
                        signif(as.data.frame(resLFCOrdered.treat),5))

write.table(out.treat, outf.treat, sep='\t', quote=F, col.names = T, row.names = F)


###############################################################################
#     genotype in untreated mice      ####################################################
###############################################################################
dds <- DESeqDataSetFromMatrix(countData=cts[, rownames(coldata)[coldata$condition=='Vehicle']],
                              colData=coldata[coldata$condition=='Vehicle',],
                              design= ~ genotype)
# Set reference factor levels
dds$genotype <- relevel(dds$genotype, ref = "Bl6_wt")
# only keep rows >= 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# DO IT
dds <- DESeq(dds, parallel=parallel)
resultsNames(dds)
# PCA PLOT
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='Mouse RNA-seq PCA plot')

# dispersion plot
plotDispEsts(dds)

test.condition.genotype <- 'genotype_JAK2V617F_vs_Bl6_wt'
res.genotype <- results(dds, name=test.condition.genotype, alpha = 0.01,)
resOrdered.genotype <- res.genotype[order(res.genotype$pvalue),]
resLFC.genotype <- lfcShrink(dds, type="ashr", res = res.genotype)
resLFCOrdered.genotype <- resLFC.genotype[order(resLFC.genotype$pvalue),]
summary(resLFC.genotype)
plotMA(resLFC.genotype, ylim=c(-3,3), main="Genotype effect in Vehicle treated")

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 9
plot.genes <- counts(dds)[order(resLFC.genotype$padj)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$genotype <- factor(coldata[pgm$sample,"genotype"], levels=c('Bl6_wt', 'JAK2V617F'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))

ggplot(pgm, aes(x=genotype, y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75)) +
    labs(title = 'Genotype effect in Vehicle treated')

outf.genotype <- file.path(outdir, 'diffex_genotype_VEHICLE_ONLY.tsv')
out.genotype <- data.frame(name = rownames(resLFCOrdered.genotype),
                           signif(as.data.frame(resLFCOrdered.genotype),5))

write.table(out.genotype, outf.genotype, sep='\t', quote=F, col.names = T, row.names = F)

###############################################################################
#     genotype in treated mice      ####################################################
###############################################################################
dds <- DESeqDataSetFromMatrix(countData=cts[, rownames(coldata)[coldata$condition=='LSDi']],
                              colData=coldata[coldata$condition=='LSDi',],
                              design= ~ genotype)
# Set reference factor levels
dds$genotype <- relevel(dds$genotype, ref = "Bl6_wt")
# only keep rows >= 10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# DO IT
dds <- DESeq(dds, parallel=parallel)
resultsNames(dds)
# PCA PLOT
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=genotype)) +
    geom_point(size=5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    theme_bw() +
    labs(title='Mouse RNA-seq PCA plot')

# dispersion plot
plotDispEsts(dds)

test.condition.genotype <- 'genotype_JAK2V617F_vs_Bl6_wt'
res.genotype <- results(dds, name=test.condition.genotype, alpha = 0.01,)
resOrdered.genotype <- res.genotype[order(res.genotype$pvalue),]
resLFC.genotype <- lfcShrink(dds, type="ashr", res = res.genotype)
resLFCOrdered.genotype <- resLFC.genotype[order(resLFC.genotype$pvalue),]
summary(resLFC.genotype)
plotMA(resLFC.genotype, ylim=c(-3,3), main="Genotype effect in LSDi treated")

# Look at the top N diffex genes
# and do a boxplot for them all in one pane
plot.n <- 9
plot.genes <- counts(dds)[order(resLFC.genotype$padj)[1:plot.n],]
pgm <- melt(plot.genes,as.is = T)
colnames(pgm) <- c('gene', 'sample', 'normalized count')
pgm$genotype <- factor(coldata[pgm$sample,"genotype"], levels=c('Bl6_wt', 'JAK2V617F'))
pgm$gene <- factor(pgm$gene, levels=rownames(plot.genes))

ggplot(pgm, aes(x=genotype, y=`normalized count`)) +
    geom_boxplot() +
    facet_wrap(.~gene, scales = 'free_y') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75)) +
    labs(title = 'Genotype effect in Vehicle treated')

outf.genotype <- file.path(outdir, 'diffex_genotype_LSDi_ONLY.tsv')
out.genotype <- data.frame(name = rownames(resLFCOrdered.genotype),
                           signif(as.data.frame(resLFCOrdered.genotype),5))

write.table(out.genotype, outf.genotype, sep='\t', quote=F, col.names = T, row.names = F)
