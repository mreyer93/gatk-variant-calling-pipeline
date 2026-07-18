# try with the ensembl gene names instead to see if I'm missing anything
library(dplyr)
library(msigdbr)
library(clusterProfiler)
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

count.f <- "~/pcloud_sync/project/rnaseq/jak2_experiment/results/counts/all.tsv"
cts <- read.table(count.f, header=TRUE, row.names="gene", check.names=FALSE)
# change rownames of count table first
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene.ids <- rownames(cts)
foo.gene <- getBM(attributes=c('ensembl_gene_id',
                               'external_gene_name'),
                  filters = 'ensembl_gene_id',
                  values = gene.ids,
                  mart = ensembl)

ens.to.gene <- setNames(foo.gene$external_gene_name, foo.gene$ensembl_gene_id)

indir <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/ensembl_genes/'
inf.treat <- file.path(indir, 'diffex_treatment.tsv')
inf.treat.jak2 <- file.path(indir, 'diffex_treatment_JAK2V617F_ONY.tsv')
inf.treat.wt <- file.path(indir, 'diffex_treatment_WT_ONY.tsv')
inf.genotype <- file.path(indir, 'diffex_genotype.tsv')
inf.genotype.vehicle <- file.path(indir, 'diffex_genotype_VEHICLE_ONLY.tsv')
inf.genotype.lsdi <- file.path(indir, 'diffex_genotype_LSDi_ONLY.tsv')
inf.interaction <- file.path(indir, 'diffex_interaction.tsv')
outdir <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/GSEA'


diffex.treat <- read.table(inf.treat, sep='\t', quote='', header = T)
diffex.treat.jak2 <- read.table(inf.treat.jak2, sep='\t', quote='', header = T)
diffex.treat.wt <- read.table(inf.treat.wt, sep='\t', quote='', header = T)
diffex.genotype <- read.table(inf.genotype, sep='\t', quote='', header = T)
diffex.genotype.vehicle <- read.table(inf.genotype.vehicle, sep='\t', quote='', header = T)
diffex.genotype.lsdi <- read.table(inf.genotype.lsdi, sep='\t', quote='', header = T)
diffex.interaction <- read.table(inf.interaction, sep='\t', quote='', header = T)
df.list <- list(diffex.treat, diffex.treat.jak2, diffex.treat.wt, diffex.genotype,diffex.genotype.vehicle,diffex.genotype.lsdi, diffex.interaction)
names(df.list) <- c('treat', 'treat_jak2', 'treat_wt', 'genotype', 'genotype_vehicle', 'genotype_lsdi', 'interaction')

H_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, ensembl_gene)
C1_t2g <- msigdbr(species = "Mus musculus", category = "C1") %>% dplyr::select(gs_name, ensembl_gene)
C2_t2g <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::select(gs_name, ensembl_gene)
C3_t2g <- msigdbr(species = "Mus musculus", category = "C3") %>% dplyr::select(gs_name, ensembl_gene)
C4_t2g <- msigdbr(species = "Mus musculus", category = "C4") %>% dplyr::select(gs_name, ensembl_gene)
C5_t2g <- msigdbr(species = "Mus musculus", category = "C5") %>% dplyr::select(gs_name, ensembl_gene)
C6_t2g <- msigdbr(species = "Mus musculus", category = "C6") %>% dplyr::select(gs_name, ensembl_gene)
C7_t2g <- msigdbr(species = "Mus musculus", category = "C7") %>% dplyr::select(gs_name, ensembl_gene)
C8_t2g <- msigdbr(species = "Mus musculus", category = "C8") %>% dplyr::select(gs_name, ensembl_gene)

# function to convert a list of ensembl genes to symbols from GSEA result
ens_str_to_symbols <- function(a){
    a.list <- strsplit(a, split='/')[[1]]
    a.sym <- ens.to.gene[a.list]
    return(paste(a.sym, collapse='/'))
}

for (n in names(df.list)){
    print(n)
    outdir.n <- file.path(outdir, n)
    dir.create(outdir.n,showWarnings = F)
    use.diffex <- df.list[[n]]
    geneList <- use.diffex$log2FoldChange
    names(geneList) <- as.character(use.diffex$name)
    geneList <- sort(geneList, decreasing = TRUE)
    gene <- names(geneList)[abs(geneList) > 2]

    # em <- enricher(gene, TERM2GENE=H_t2g, pvalueCutoff = 0.1)
    em_H <- GSEA(geneList, TERM2GENE = H_t2g, pvalueCutoff = 0.1)
    em_C1 <- GSEA(geneList, TERM2GENE = C1_t2g, pvalueCutoff = 0.1)
    em_C2 <- GSEA(geneList, TERM2GENE = C2_t2g, pvalueCutoff = 0.1)
    em_C3 <- GSEA(geneList, TERM2GENE = C3_t2g, pvalueCutoff = 0.1)
    em_C4 <- GSEA(geneList, TERM2GENE = C4_t2g, pvalueCutoff = 0.1)
    em_C5 <- GSEA(geneList, TERM2GENE = C5_t2g, pvalueCutoff = 0.1)
    em_C6 <- GSEA(geneList, TERM2GENE = C6_t2g, pvalueCutoff = 0.1)
    em_C7 <- GSEA(geneList, TERM2GENE = C7_t2g, pvalueCutoff = 0.1)
    em_C8 <- GSEA(geneList, TERM2GENE = C8_t2g, pvalueCutoff = 0.1)

    # Core enriched genes: This is the subset of genes that contributes most to the enrichment result
    em_H@result$core_enrichment_symbol <- sapply(em_H@result$core_enrichment, function(x) ens_str_to_symbols(x))
    em_C1@result$core_enrichment_symbol <- sapply(em_C1@result$core_enrichment, function(x) ens_str_to_symbols(x))
    em_C2@result$core_enrichment_symbol <- sapply(em_C2@result$core_enrichment, function(x) ens_str_to_symbols(x))
    em_C3@result$core_enrichment_symbol <- sapply(em_C3@result$core_enrichment, function(x) ens_str_to_symbols(x))
    em_C4@result$core_enrichment_symbol <- sapply(em_C4@result$core_enrichment, function(x) ens_str_to_symbols(x))
    em_C5@result$core_enrichment_symbol <- sapply(em_C5@result$core_enrichment, function(x) ens_str_to_symbols(x))
    em_C6@result$core_enrichment_symbol <- sapply(em_C6@result$core_enrichment, function(x) ens_str_to_symbols(x))
    em_C7@result$core_enrichment_symbol <- sapply(em_C7@result$core_enrichment, function(x) ens_str_to_symbols(x))
    em_C8@result$core_enrichment_symbol <- sapply(em_C8@result$core_enrichment, function(x) ens_str_to_symbols(x))

    # em.df.list <- list(em_H @result[,1:min(ncol(em_H@result), 10)], em_C1@result[,1:min(ncol(em_C1@result), 10)], em_C2@result[,1:min(ncol(em_C2@result), 10)], em_C3@result[,1:min(ncol(em_C3@result), 10)], em_C4@result[,1:min(ncol(em_C4@result), 10)], em_C5@result[,1:min(ncol(em_C5@result), 10)], em_C6@result[,1:min(ncol(em_C6@result), 10)], em_C7@result[,1:min(ncol(em_C7@result), 10)], em_C8@result[,1:min(ncol(em_C8@result), 10)])
    em.df.list <- list(em_H@result, em_C1@result, em_C2@result, em_C3@result, em_C4@result, em_C5@result, em_C6@result, em_C7@result, em_C8@result)
    names(em.df.list) <- c("msidgdb_H", "msidgdb_C1", "msidgdb_C2", "msidgdb_C3", "msidgdb_C4", "msidgdb_C5", "msidgdb_C6", "msidgdb_C7", "msidgdb_C8")

    for (m in names(em.df.list)){
        outf <- file.path(outdir.n, paste0(m, '.tsv'))
        write.table(em.df.list[[m]], outf, sep='\t', quote = F, row.names = F, col.names = T)
    }

}


# plot(gvisTable(em_H@result[,1:10]))
# # plot(gvisTable(em_C1@result[,1:10]))
# plot(gvisTable(em_C2@result[,1:10]))
# plot(gvisTable(em_C3@result[,1:10]))
# plot(gvisTable(em_C4@result[,1:10]))
# plot(gvisTable(em_C5@result[,1:10]))
# plot(gvisTable(em_C6@result[,1:10]))
# plot(gvisTable(em_C7@result[,1:10]))
# plot(gvisTable(em_C8@result[,1:10]))
#
# # H   hallmark gene sets  are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.
# # C1  positional gene sets  for each human chromosome and cytogenetic band.
# # C2  curated gene sets  from online pathway databases, publications in PubMed, and knowledge of domain experts.
# # C3  regulatory target gene sets  based on gene target predictions for microRNA seed sequences and predicted transcription factor binding sites.
# # C4  computational gene sets  defined by mining large collections of cancer-oriented microarray data.
# # C5  ontology gene sets  consist of genes annotated by the same ontology term.
# # C6  oncogenic signature gene sets  defined directly from microarray gene expression data from cancer gene perturbations.
# # C7  immunologic signature gene sets  represent cell states and perturbations within the immune system.
# # C8  cell type signature gene sets  curated from cluster markers identified in single-cell sequencing studies of human tissue.


gseaplot(em_H, 'HALLMARK_IL6_JAK_STAT3_SIGNALING', title = 'HALLMARK_IL6_JAK_STAT3_SIGNALING, JAK2 mutant in untreated')


# 2021-06-16
# Make RNA-Seq figure for Ria
# bar graph with NES score, for C3
#  - I want to highlight PU1/ETS
#  - It is the top hit even though PU1 itself is not downregulated or upregulated by LSD1i treatment
#  - WT cells (comparing treated and untreated) only have MALM1 as hit for your C3 analysis
# THIS IS IN MUTANT CELLS ONLY

res <- em_C3@result %>% arrange(p.adjust, desc(abs(NES)))
res <- res[]

n <- 12
plotlist <- lapply(1:n, function(i) gseaplot(em_C3, res[i,"ID"], cex_label_category=1, title =  paste0(res[i,"ID"], '   NES=', round(res[i,"NES"],2))))

pdf('~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/GSEA/treat_jak2/treat_jak2_C3_top.pdf', height=16, width=14, pointsize = 6)
cowplot::plot_grid(plotlist=plotlist[1:6], nrow=3)
cowplot::plot_grid(plotlist=plotlist[7:12], nrow=3)
dev.off()

# just a barplot?
rp <- res[1:12,]
rp$ID <- factor(rp$ID, levels=rev(rp$ID))

pdf('~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/GSEA/treat_jak2/treat_jak2_C3_barplot.pdf', height=6, width=5)
ggplot(rp, aes(x=NES, y=ID)) +
    geom_bar(stat='identity', fill='steelblue') +
    theme_bw() +
    labs(title='Top differential C3 gene sets',
         subtitle='JAK2 mutant, LSDi mice',
         x='Normalized enrichment score',
         y='Gene Set')
dev.off()

