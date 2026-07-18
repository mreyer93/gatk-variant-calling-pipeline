# get the gene set information for Ria that she asked for

# Can you pull for me the log-fold change and q-values
# for the genes included in the different gene sets
# (HALLMARK INFA and INFG) for
# JAK2VF_Genotype_VEHICLE ONLY, JAK2_TREATMENT
# MPL_Genotype, MPL_TREATMENT
# (should be about 200 for IFNG and 176 for IFNA)
# I want to have side by side comparison for GENOTYPE and TREATMENT.

library(msigdbr)
library(clusterProfiler)
library(biomaRt)
library(dplyr)


# read in differential expression analysis
inf.jak2.genotype.veh <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_genotype_VEHICLE_ONLY.tsv'
inf.jak2.treat <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/diffex_treatment.tsv'
diffex.jak2.genotype.veh <- read.table(inf.jak2.genotype.veh, sep='\t', quote='', header = T)
diffex.jak2.treat <- read.table(inf.jak2.treat, sep='\t', quote='', header = T)
inf.mpl.genotype <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex/diffex_genotype_SRA.tsv'
inf.mpl.treat <- '~/pcloud_sync/project/rnaseq/mpl_experiment/results/diffex/diffex_treatment_BOX.tsv'
diffex.mpl.genotype <- read.table(inf.mpl.genotype, sep='\t', quote='', header = T)
diffex.mpl.treat <- read.table(inf.mpl.treat, sep='\t', quote='', header = T)

rownames(diffex.jak2.genotype.veh) <- diffex.jak2.genotype.veh$name
rownames(diffex.jak2.treat) <- diffex.jak2.treat$name
rownames(diffex.mpl.genotype) <- diffex.mpl.genotype$name
rownames(diffex.mpl.treat) <- diffex.mpl.treat$name


# rowname map
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene.ids <- unique(c(diffex.jak2.genotype.veh$name, diffex.jak2.treat$name))
foo.gene <- getBM(attributes=c('ensembl_gene_id',
                               'external_gene_name'),
                  filters = 'ensembl_gene_id',
                  values = gene.ids,
                  mart = ensembl)
ens.to.gene <- setNames(foo.gene$external_gene_name, foo.gene$ensembl_gene_id)




H_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% dplyr::select(gs_name, gene_symbol)

gs1.name <- 'HALLMARK_INTERFERON_ALPHA_RESPONSE'
gs2.name <- 'HALLMARK_INTERFERON_GAMMA_RESPONSE'
gs1.symbol <- pull(H_t2g[H_t2g$gs_name==gs1.name,], gene_symbol)
gs2.symbol <- pull(H_t2g[H_t2g$gs_name==gs2.name,"gene_symbol"], gene_symbol)

df.alpha <- data.frame(gene=gs1.symbol,
                       gene.set = gs1.name,
                       lfc.jak2.genotype.veh = diffex.jak2.genotype.veh[gs1.symbol, "log2FoldChange"],
                       qv.jak2.genotype.veh = diffex.jak2.genotype.veh[gs1.symbol, "padj"],
                       lfc.jak2.treat = diffex.jak2.treat[gs1.symbol, "log2FoldChange"],
                       qv.jak2.treat = diffex.jak2.treat[gs1.symbol, "padj"],
                       lfc.mpl.genotype = diffex.mpl.genotype[gs1.symbol, "log2FoldChange"],
                       qv.mpl.genotype = diffex.mpl.genotype[gs1.symbol, "padj"],
                       lfc.mpl.treat = diffex.mpl.treat[gs1.symbol, "log2FoldChange"],
                       qv.mpl.treat = diffex.mpl.treat[gs1.symbol, "padj"]
                       )

df.gamma <- data.frame(gene=gs2.symbol,
                       gene.set = gs2.name,
                       lfc.jak2.genotype.veh = diffex.jak2.genotype.veh[gs2.symbol, "log2FoldChange"],
                       qv.jak2.genotype.veh = diffex.jak2.genotype.veh[gs2.symbol, "padj"],
                       lfc.jak2.treat = diffex.jak2.treat[gs2.symbol, "log2FoldChange"],
                       qv.jak2.treat = diffex.jak2.treat[gs2.symbol, "padj"],
                       lfc.mpl.genotype = diffex.mpl.genotype[gs2.symbol, "log2FoldChange"],
                       qv.mpl.genotype = diffex.mpl.genotype[gs2.symbol, "padj"],
                       lfc.mpl.treat = diffex.mpl.treat[gs2.symbol, "log2FoldChange"],
                       qv.mpl.treat = diffex.mpl.treat[gs2.symbol, "padj"]
                       )

write.table(df.alpha, '~/pcloud_sync/project/rnaseq/interferon_alpha_diffex.tsv', sep='\t', quote=F, row.names = F, col.names = T)
write.table(df.gamma, '~/pcloud_sync/project/rnaseq/interferon_gamma_diffex.tsv', sep='\t', quote=F, row.names = F, col.names = T)
