library(msigdbr)
library(clusterProfiler)

# read in differential expression analysis
indir <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex'
inf.treat <- file.path(indir, 'diffex_treatment.tsv')
inf.treat.jak2 <- file.path(indir, 'diffex_treatment_JAK2V617F_ONY.tsv')
inf.treat.wt <- file.path(indir, 'diffex_treatment_WT_ONY.tsv')
inf.genotype <- file.path(indir, 'diffex_genotype.tsv')
inf.interaction <- file.path(indir, 'diffex_interaction.tsv')

diffex.treat <- read.table(inf.treat, sep='\t', quote='', header = T)
diffex.treat.jak2 <- read.table(inf.treat.jak2, sep='\t', quote='', header = T)
diffex.treat.wt <- read.table(inf.treat.wt, sep='\t', quote='', header = T)
diffex.genotype <- read.table(inf.genotype, sep='\t', quote='', header = T)
diffex.interaction <- read.table(inf.interaction, sep='\t', quote='', header = T)

# all_gene_sets = msigdbr(species = "Mus musculus") %>%
#     dplyr::select(gs_name, gene_symbol)
# head(all_gene_sets)
# msigdbr_collections()


H_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, gene_symbol)
C1_t2g <- msigdbr(species = "Mus musculus", category = "C1") %>% dplyr::select(gs_name, gene_symbol)
C2_t2g <- msigdbr(species = "Mus musculus", category = "C2") %>% dplyr::select(gs_name, gene_symbol)
C3_t2g <- msigdbr(species = "Mus musculus", category = "C3") %>% dplyr::select(gs_name, gene_symbol)
C4_t2g <- msigdbr(species = "Mus musculus", category = "C4") %>% dplyr::select(gs_name, gene_symbol)
C5_t2g <- msigdbr(species = "Mus musculus", category = "C5") %>% dplyr::select(gs_name, gene_symbol)
C6_t2g <- msigdbr(species = "Mus musculus", category = "C6") %>% dplyr::select(gs_name, gene_symbol)
C7_t2g <- msigdbr(species = "Mus musculus", category = "C7") %>% dplyr::select(gs_name, gene_symbol)
C8_t2g <- msigdbr(species = "Mus musculus", category = "C8") %>% dplyr::select(gs_name, gene_symbol)

use.diffex <- diffex.treat
geneList <- use.diffex$log2FoldChange

## feature 2: named vector
names(geneList) <- as.character(use.diffex$name)
# feature 3: decreasing order
geneList <- sort(geneList, decreasing = TRUE)

gene <- names(geneList)[abs(geneList) > 2]
head(gene)

# em <- enricher(gene, TERM2GENE=m_t2g, pvalueCutoff = 1)
em_H <- GSEA(geneList, TERM2GENE = H_t2g, pvalueCutoff = 1)
em_C1 <- GSEA(geneList, TERM2GENE = C1_t2g, pvalueCutoff = 1)
em_C2 <- GSEA(geneList, TERM2GENE = C2_t2g, pvalueCutoff = 1)
em_C3 <- GSEA(geneList, TERM2GENE = C3_t2g, pvalueCutoff = 1)
em_C4 <- GSEA(geneList, TERM2GENE = C4_t2g, pvalueCutoff = 1)
em_C5 <- GSEA(geneList, TERM2GENE = C5_t2g, pvalueCutoff = 1)
em_C6 <- GSEA(geneList, TERM2GENE = C6_t2g, pvalueCutoff = 1)
em_C7 <- GSEA(geneList, TERM2GENE = C7_t2g, pvalueCutoff = 1)
em_C8 <- GSEA(geneList, TERM2GENE = C8_t2g, pvalueCutoff = 1)

head(em_H)
head(em_C1)
head(em_C2)
head(em_C3)
head(em_C4)
head(em_C5)
head(em_C6)
head(em_C7)
head(em_C8)
c(min(em_H$p.adjust), min(em_C1$p.adjust), min(em_C2$p.adjust), min(em_C3$p.adjust), min(em_C4$p.adjust), min(em_C5$p.adjust), min(em_C6$p.adjust), min(em_C7$p.adjust), min(em_C8$p.adjust))
