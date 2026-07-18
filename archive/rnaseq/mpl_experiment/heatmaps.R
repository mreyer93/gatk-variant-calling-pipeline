# Making heatmaps for the RNA-Seq data
# in the MPL mutant, IMG7289
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(viridis)
library(pals)

count.f <- "~/pcloud_sync/project/rnaseq/mpl_experiment/results/counts/all.tsv"
norm.count.sra.f <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/counts/normalized_counts_export_SRA.tsv'
norm.count.box.f <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/counts/normalized_counts_export_BOX.tsv'

sample.f <- "~/pcloud_sync/project/rnaseq/mpl_experiment/config/samples.tsv"
units.f <- "~/pcloud_sync/project/rnaseq/mpl_experiment/config/units.tsv"
outdir <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex'

# metadata
coldata <- read.table(sample.f, header=TRUE, row.names="sample_name", check.names=FALSE)
coldata.sra <- coldata[coldata$batch == 'SRA',]
coldata.box <- coldata[coldata$batch == 'BOX',]

# change sample order
if (F){
sample.order <- c("9292_wt_vehicle_MEP","9293_wt_vehicle_MEP","9373_wt_vehicle_MEP",
                  "JAF8_J2VF_vehicle_MEP","JAG1_J2VF_vehicle_MEP","JAG7_J2VF_vehicle_MEP",
                  "9198_wt_LSDi_MEP","9374_wt_LSDi_MEP","JAJ3_wt_LSDi_MEP",
                  "JAE3_J2VF_LSDi_MEP","JAF9_J2VF_LSDi_MEP","JAG3_J2VF_LSDi_MEP")
coldata <- coldata[sample.order, ]
}

# raw and normalized counts
cts <- read.table(count.f, header=TRUE, row.names="gene", check.names=FALSE)
counts.norm.sra <- read.table(norm.count.sra.f, sep='\t', quote='', header=T, check.names = F)
counts.norm.sra.mat <- counts.norm.sra[,3:ncol(counts.norm.sra)]
counts.norm.box <- read.table(norm.count.box.f, sep='\t', quote='', header=T, check.names = F)
counts.norm.box.mat <- counts.norm.box[,3:ncol(counts.norm.box)]

rownames(counts.norm.sra.mat) <- counts.norm.sra$name
rownames(counts.norm.box.mat) <- counts.norm.box$name

# DEGs
input.df <- data.frame(files= c("/path/to/data/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex/diffex_genotype_SRA.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex/diffex_treatment_BOX.tsv"
                                ),
                       names = c("genotype_SRA",
                                 "treatment_BOX"
                                 )
)

pdf('~/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex/heatmaps/combined.pdf', height=10, width=8)
cluster.outf <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex/heatmaps/cluster_assignments.tsv'
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

    if (name=='genotype_SRA'){
        cnm.plot <- counts.norm.sra.mat[dex.gene.list, ]
        rownames(cnm.plot.means) <- rownames(cnm.plot)
    } else if (name=='treatment_BOX'){
        cnm.plot <- counts.norm.box.mat[dex.gene.list, ]
    }
    cnm.plot.scale <- t(scale(t(cnm.plot)))

    # aggregate each of these to mean to add to the gene table
    cnm.plot.means <- data.frame(group_1_mean=rowMeans(cnm.plot[,1:3]),
                                 group_2_mean=rowMeans(cnm.plot[,4:6]))

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
