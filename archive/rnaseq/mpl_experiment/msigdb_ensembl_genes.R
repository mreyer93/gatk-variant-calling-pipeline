# MSIGDB profiles for the MPLW515L experiment
library(dplyr)
library(msigdbr)
library(biomaRt)
library(clusterProfiler)


indir <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex/ensembl_genes/'
inf.treat <- file.path(indir, 'diffex_treatment_BOX.tsv')
inf.genotype <- file.path(indir, 'diffex_genotype_SRA.tsv')
outdir <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex/GSEA'


diffex.treat <- read.table(inf.treat, sep='\t', quote='', header = T)
diffex.genotype <- read.table(inf.genotype, sep='\t', quote='', header = T)
df.list <- list(diffex.treat, diffex.genotype)
names(df.list) <- c('treat', 'genotype')

# rowname map
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene.ids <- unique(c(diffex.treat$name, diffex.genotype$name))
foo.gene <- getBM(attributes=c('ensembl_gene_id',
                               'external_gene_name'),
                  filters = 'ensembl_gene_id',
                  values = gene.ids,
                  mart = ensembl)
ens.to.gene <- setNames(foo.gene$external_gene_name, foo.gene$ensembl_gene_id)

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
    dir.create(outdir.n, recursive = T, showWarnings = F)
    use.diffex <- df.list[[n]]
    geneList <- use.diffex$log2FoldChange
    names(geneList) <- as.character(use.diffex$name)
    geneList <- sort(geneList, decreasing = TRUE)
    gene <- names(geneList)[abs(geneList) > 2]

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


library(cowplot)
# treatment sets
plot.sets <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION" ,"HALLMARK_UV_RESPONSE_DN" ,"HALLMARK_DNA_REPAIR" ,"HALLMARK_HEME_METABOLISM")
plotlist  <- lapply(plot.sets, function(x) gseaplot(em_H, x, title = paste0(x, '\n IMG7289 treatment in MPLW515L mice')))
plot_grid(plotlist = plotlist, ncol=2)

# genotype sets
plot.sets <- c("HALLMARK_P53_PATHWAY" ,"HALLMARK_ESTROGEN_RESPONSE_EARLY" ,"HALLMARK_MTORC1_SIGNALING" ,"HALLMARK_TNFA_SIGNALING_VIA_NFKB" )
plotlist  <- lapply(plot.sets, function(x) gseaplot(em_H, x, title = paste0(x, '\n MPLW515L mutation in untreated mice')))
plot_grid(plotlist = plotlist, ncol=2)
