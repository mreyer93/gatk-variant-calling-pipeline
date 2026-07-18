library(clusterProfiler)
library(biomaRt)
library(pathview)
httr::set_config(httr::config(ssl_verifypeer = FALSE))

# investigation of KEGG pathways with the differential expression data
# DEGs
# this will be looking at genes that are changed by LSDi in JAK2 mutant background
input.df <- data.frame(files= c("/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/ensembl_genes/diffex_genotype_LSDi_ONLY.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/ensembl_genes/diffex_genotype.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/ensembl_genes/diffex_genotype_VEHICLE_ONLY.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/ensembl_genes/diffex_interaction.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/ensembl_genes/diffex_treatment_JAK2V617F_ONY.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/ensembl_genes/diffex_treatment.tsv",
                                "/path/to/data/pcloud_sync/project/rnaseq/jak2_experiment/results/diffex/ensembl_genes/diffex_treatment_WT_ONY.tsv"
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
outdir.tables <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/kegg_pathways/tables'
outdir.images <- '~/pcloud_sync/project/rnaseq/jak2_experiment/results/kegg_pathways/images'

diffex.list <- lapply(input.df$files, function(f) read.table(f, header=T, sep='\t', quote=''))
names(diffex.list) <- input.df$name
all.gene.names <- unique(unlist(lapply(diffex.list, function(x) x$name)))

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

# convert and add entrez gene id
foo.gene <- getBM(attributes=c('ensembl_gene_id',
                               'entrezgene_id'),
                  values = all.gene.names,
                  mart = ensembl)
foo.geneName <- getBM(attributes=c('entrezgene_id',
                               'external_gene_name'),
                  values = all.gene.names,
                  mart = ensembl)

ensembl.to.entrez <- setNames(foo.gene$entrezgene_id, foo.gene$ensembl_gene_id)
entrez.to.name <- setNames(foo.geneName$external_gene_name,foo.geneName$entrezgene_id)

diffex.list2 <- lapply(diffex.list, function(x) {
    x$entrez <- ensembl.to.entrez[x$name]
    print(sum(is.na(x$entrez)))
    return(x[!is.na(x$entrez),])
})


# do for each comparison
for (comparison in input.df$names){
    print('STARTING #######################')
    print(comparison)
    padj.cutoff <- 0.01
    abs.lfc.cutoff <- 1
    sig.genes <- diffex.list2[[comparison]]$entrez[diffex.list2[[comparison]]$padj < padj.cutoff &
                                              abs(diffex.list2[[comparison]]$log2FoldChange) > abs.lfc.cutoff]

    kk <- enrichKEGG(gene = sig.genes, organism = 'mmu')
    head(kk, n=10)
    kkdf <- as.data.frame(kk)
    # add gene names back in
    kkdf$gene_names <- sapply(kkdf$geneID, function(x) paste(entrez.to.name[strsplit(x, '/')[[1]]], collapse = '/'))

    # write out a table for this
    outf.table <- file.path(outdir.tables, paste(comparison, '.tsv'))
    write.table(kkdf, outf.table, sep='\t', quote=F, row.names = F, col.names = T)

    outdir.images.comparison <- file.path(outdir.images, comparison)
    dir.create(outdir.images.comparison, recursive = T, showWarnings = F)

    # make figures on the pathways
    setwd(outdir.images.comparison)
    # only highlight images that are significant
    this.df <- diffex.list2[[comparison]][diffex.list2[[comparison]]$padj <= padj.cutoff,]
    logFC <- this.df$log2FoldChange
    names(logFC) <- this.df$entrez

    for (pick.pathway in kkdf$ID){
        print(pick.pathway)
        pathway.name <- gsub(" ", "_", kkdf[pick.pathway, "Description"])
        pathview(gene.data = logFC,
                 pathway.id = pick.pathway,
                 species = "mmu",
                 kegg.dir = '~/pcloud_sync/project/rnaseq/jak2_experiment/results/kegg_pathways/images/',
                 out.suffix = paste(pathway.name, comparison, sep='__'),
                 limit = list(gene=3, cpd=1),
                 low= list(gene='steelblue', cpd='steelblue'))

    }
}

# can I make a heatmap of the top pathways per experimental condition?
# read the tables back in...
top.df <- do.call(rbind, lapply(input.df$names, function(comparison) {
    inf <- file.path(outdir.tables, paste(comparison, '.tsv'))
    df <- read.table(inf, sep='\t', quote='', header=T)
    df$comparison <- comparison
    return(df)
}))
