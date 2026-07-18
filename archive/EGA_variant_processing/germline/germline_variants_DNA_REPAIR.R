# association of germline variants with LOH, in the EGA data
library(ggplot2)
library(viridis)
library(gplots)
library(reshape2)
library(msigdbr)
library(ggrepel)
library(ggpp)
vcf.20 <- read.table('~/pcloud_sync/project/EGA_data/results/germline/vcf_20.tsv', sep='\t', quote='', header = T)
vcf.20.meta <- read.table('~/pcloud_sync/project/EGA_data/results/germline/vcf_20_meta.tsv', sep='\t', quote='', header = T)
# top 100 different variants, high in EGA but low in gnomAD
vcf.20.topdiff <- vcf.20[vcf.20$gnomad.diff >= 0.05, ]
vcf.20.topdiff <- vcf.20.topdiff[order(vcf.20.topdiff$gnomad.l2fc, decreasing = T)[1:min(1000, nrow(vcf.20.topdiff))], ]

remove.genes <- c('HYDIN')
# vcf.20.topdiff <- vcf.20.topdiff[,!(colnames(vcf.20.topdiff) %in% remove.cols)]
vcf.20.topdiff <- vcf.20.topdiff[!(vcf.20.topdiff$hugoSymbol %in% remove.genes),]
# eliminate things found in >90% of patients because they're surely false positives
vcf.20.topdiff <- vcf.20.topdiff[vcf.20.topdiff$frac.pt.with.variant < 0.75, ]
# less thatn 20% in gnomad
vcf.20.topdiff <- vcf.20.topdiff[vcf.20.topdiff$gnomad.to.plot < 0.35, ]



H_t2g <- as.data.frame(msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, human_gene_symbol))
dna.name <- 'HALLMARK_DNA_REPAIR'
dna.genes <- H_t2g[H_t2g$gs_name==dna.name, "human_gene_symbol"]

vcf.20.topdiff[vcf.20.topdiff$hugoSymbol %in% dna.genes, ]
vcf.20.dna <- vcf.20[vcf.20$hugoSymbol %in% dna.genes, ]

ggplot(vcf.20.dna, aes(x=gnomad.to.plot, y=frac.alleles, label=gene_change)) +
    geom_point() +
    theme_bw() +
    labs(title='Variants in DNA repair genes',
         x='gnomAD fraction of alleles', y='EGA fraction of alleles') +
    stat_dens2d_labels(geom = "text_repel", keep.fraction = 0.25)


# what about recombination gene set from Hugh?
# Added all the DNA polymerases
rec.genes <- c('RAD50', 'RAD51', 'MRE1', 'NBS1', 'ATM', 'RECA',
               "DNTT", "POLA1", "POLA2", "POLB", "POLD1", "POLD2", "POLD3", "POLD4", "POLE", "POLE2", "POLE3", "POLE4", "POLG", "POLG2", "POLH", "POLI", "POLK", "POLL", "POLM", "POLN", "POLQ", "REV3L", "MAD2L2", "PRIMPOL", "REV1",
               'BRCA1', 'BRCA2')

vcf.20.rec <- vcf.20[vcf.20$hugoSymbol %in% rec.genes, ]
vcf.20.rec <- vcf.20.rec[order(vcf.20.rec$gnomad.l2fc, decreasing = T), ]

ggplot(vcf.20.rec, aes(x=gnomad.to.plot, y=frac.alleles, label=gene_change)) +
    geom_point() +
    theme_bw() +
    labs(title='Variants in recombination genes',
         x='gnomAD fraction of alleles', y='EGA fraction of alleles') +
    stat_dens2d_labels(geom = "text_repel", keep.fraction = 0.25)  +
    ylim(c(NA, 1)) +
    xlim(c(NA, 1))
