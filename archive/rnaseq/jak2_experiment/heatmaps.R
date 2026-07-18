# Making heatmaps for the RNA-Seq data
# starting with batch 1
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(viridis)
library(pals)

count.f <- "~/pcloud_sync/project/rnaseq/jak2_experiment/results/counts/all.tsv"
norm.count.f <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/counts/normalized_counts_export.tsv'

sample.f <- "~/pcloud_sync/project/rnaseq/jak2_experiment/config/samples.tsv"
units.f <- "~/pcloud_sync/project/rnaseq/jak2_experiment/config/units.tsv"
outdir <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex'

# metadata
coldata <- read.table(sample.f, header=TRUE, row.names="sample_name", check.names=FALSE)
# change sample order
sample.order <- c("9292_wt_vehicle_MEP","9293_wt_vehicle_MEP","9373_wt_vehicle_MEP",
                  "JAF8_J2VF_vehicle_MEP","JAG1_J2VF_vehicle_MEP","JAG7_J2VF_vehicle_MEP",
                  "9198_wt_LSDi_MEP","9374_wt_LSDi_MEP","JAJ3_wt_LSDi_MEP",
                  "JAE3_J2VF_LSDi_MEP","JAF9_J2VF_LSDi_MEP","JAG3_J2VF_LSDi_MEP")
coldata <- coldata[sample.order, ]

# raw and normalized counts
cts <- read.table(count.f, header=TRUE, row.names="gene", check.names=FALSE)
counts.norm <- read.table(norm.count.f, sep='\t', quote='', header=T, check.names = F)
counts.norm.mat <- counts.norm[,3:ncol(counts.norm)]
rownames(counts.norm.mat) <- counts.norm$name
cts <- cts[,sample.order]
counts.norm.mat <- counts.norm.mat[,make.names(sample.order)]

# DEGs
# this will be looking at genes that are changed by LSDi in JAK2 mutant background
input.df <- data.frame(files= c("/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_genotype_LSDi_ONLY.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_genotype.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_genotype_VEHICLE_ONLY.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_interaction.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_treatment_JAK2V617F_ONY.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_treatment.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_treatment_WT_ONY.tsv"
                                ),
                       names = c("genotype_LSDi_ONLY",
                                 "genotype_all_samples",
                                 "genotype_VEHICLE_ONLY",
                                 "interaction_all_samples",
                                 "treatment_JAK2V617F_ONY",
                                 "treatment_all_samples",
                                 "treatment_WT_ONY"
                                 )
)

pdf('~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/heatmaps/combined.pdf', height=10, width=8)
cluster.outf <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/heatmaps/cluster_assignments.tsv'
for (i in 1:nrow(input.df)){
    dex.f <- input.df[i, 'files']
    name <- input.df[i,"names"]
    print(name)
    dex <- read.table(dex.f, sep='\t', quote='', header=T)
    dex[is.na(dex$padj), "padj"] <- 1
    rownames(dex) <- dex$name
    # pick genes that are changed with an abs value and pv cutoff
    min.abs.lfc <- 1.5
    min.p.adj <- 0.01
    dex.gene.list <- dex[abs(dex$log2FoldChange) >= min.abs.lfc & dex$padj < min.p.adj, "name"]
    print(length(dex.gene.list))

    cnm.plot <- counts.norm.mat[dex.gene.list, ]
    cnm.plot.scale <- t(scale(t(cnm.plot)))

    # aggregate each of these to mean to add to the gene table
    cnm.plot.means <- data.frame(wt_vehicle_mean=rowMeans(cnm.plot[,1:3]),
                                 JAK2_vehicle_mean=rowMeans(cnm.plot[,4:6]),
                                 wt_LSDi_mean=rowMeans(cnm.plot[,7:9]),
                                 JAK2_LSDi_mean=rowMeans(cnm.plot[,10:12]))
    rownames(cnm.plot.means) <- rownames(cnm.plot)

    # do some clustering of the genes, and get clusters
    gene.hclust <- hclust(dist(cnm.plot.scale), method='ward.D2')
    k <- 4
    clusters <- cutree(gene.hclust, k=k)
    pal <- brewer.pal(k, "Set1")
    cluster.cols <- pal[clusters]
    colmap = c("red", "blue", "green", "purple")
    cluster.df <- data.frame(name=name,
                             gene=names(clusters),
                             cluster=clusters,
                             color=colmap[clusters],
                             log2FoldChange=round(dex[names(clusters), "log2FoldChange"],3),
                             padj=signif(dex[names(clusters), "padj"],4))
    cluster.df <- cbind(cluster.df, cnm.plot.means[cluster.df$gene,])
    write.table(cluster.df, cluster.outf, sep='\t', quote=F, col.names = T, row.names = F, append=(i!=1))

    heatmap.2(cnm.plot.scale, Rowv=T, Colv=NA, trace='none', colsep=c(3,6,9),
          col=ocean.balance(64),
          margins=c(13,5),
          lhei = c(1,5),
          lwid = c(2,9),
          hclustfun = function(x) hclust(x, method = 'ward.D2'),
          key.xlab = 'Z-score', key.title = 'Color Key', labRow = NA,
          main=paste0("Top genes: ",name, "\n LFC>1.5 & FDR<0.01"),
          RowSideColors = cluster.cols)

}
dev.off()
