library(ggpubr)
library(rafalib)
library(ggrepel)
library(ggpp)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(viridis)

# combine germline vcf calls from different chromosomes
indir <- '~/local_data/EGA_data/06_GDB/annotation_combined'
chr.nums <- 1:22

vcf.orig <- do.call(rbind, lapply(chr.nums, function(i) {
    print(i)
    inf <- file.path(indir, paste0('chr', i,'_combined_filtered.vcf'))
    vcf.chr <- read.table(inf, sep='\t', header=T, quote='', check.names = F)
    return(vcf.chr)
}))
vcf <- vcf.orig
# set rownames to a big key
rownames(vcf) <- paste(vcf$CHROM, vcf$POS, vcf$REF, vcf$ALT, vcf$hugoSymbol, vcf$proteinChange, vcf$CADD_phred, sep='_')
vcf <- vcf[!is.na(vcf$CADD_phred), ]
# vcf$CADD_phred[is.na(vcf$CADD_phred)] <- 0

# some stats about the variants
nrow(vcf)
mypar(1,1,mar = c(3.5, 3.5, 1.6, 1.1))
barplot(table(vcf$CHROM)[paste0('chr', 1:22)], las=2, main="Variants per chromosome",
        ylim=c(0,50000))
hist(vcf$CADD_phred, main='CADD score distribution', xlab='CADD score')

# freq vs cadd
vcf.sample.cols <- 19:ncol(vcf)
vcf.sample.names <- colnames(vcf)[vcf.sample.cols]
# remove NA (Set to zero) in the vcf calls
vcf[,vcf.sample.names][is.na(vcf[,vcf.sample.names])] <- 0

n.pts <- length(vcf.sample.cols)
vcf.data <- vcf[,vcf.sample.cols]
# number of patients with the variant
n.pt.with.variant <- apply(vcf.data, 1, function(x) sum(x != 0, na.rm = T))
n.alleles <- apply(vcf.data, 1, function(x) sum(x, na.rm = T))
frac.pt.with.variant <- n.pt.with.variant / n.pts
hist(frac.pt.with.variant, main='Fraction of patients with variant',
     xlab='Fraction of patients')

vcf$n.pt.with.variant <- n.pt.with.variant
vcf$frac.pt.with.variant <- frac.pt.with.variant
vcf$frac.alleles <- n.alleles / (n.pts * 2)

p1 <- ggplot(vcf[vcf$frac.pt.with.variant>0,], aes(y=frac.pt.with.variant, x=CADD_phred)) +
    geom_hex(bins=70) +
    theme_bw() +
    scale_fill_continuous(type = "viridis", trans = "log10")
p2 <- ggplot(vcf[vcf$frac.pt.with.variant>0 & vcf$CADD_phred >15,], aes(y=frac.pt.with.variant, x=CADD_phred)) +
    geom_hex(bins=70) +
    theme_bw() +
    # scale_fill_continuous(type = "viridis")
    scale_fill_continuous(type = "viridis", trans = "log10")
plot_grid(p1,p2)


# cadd 20 and at least 2 patients
vcf.20 <- vcf[vcf$CADD_phred >=20 & vcf$n.pt.with.variant >=2,]
length(unique(vcf.20$hugoSymbol))

# how can I aggregate this more
vcf.20.data <- vcf.20[,vcf.sample.cols]
table(as.matrix(vcf.20.data), useNA = 'ifany')
# binary - at least heterozygotes
vcf.20.data.bin <- vcf.20.data > 0

# unique variants per gene
uniq.variants.per.gene <- aggregate(vcf.20$genomeChange, by=list(vcf.20$hugoSymbol), FUN= function(x) length(unique(x)))
uniq.variants.per.gene <- uniq.variants.per.gene[order(uniq.variants.per.gene$x, decreasing = T), ]

# incorporate the annotations from dbSNP for number for fraction of mutations in the general population
# individual variant level
# gnomad to plot, take exome, but if exome is zero, take genome
vcf.20$gnomAD_exome_AF[is.na(vcf.20$gnomAD_exome_AF)] <- 0
vcf.20$gnomAD_genome_AF[is.na(vcf.20$gnomAD_genome_AF)] <- 0
vcf.20$gnomad.to.plot <- vcf.20$gnomAD_exome_AF
vcf.20$gnomad.to.plot[vcf.20$gnomad.to.plot == 0] <- vcf.20$gnomAD_genome_AF[vcf.20$gnomad.to.plot == 0]

ggplot(vcf.20, aes(x=gnomad.to.plot, y=frac.alleles)) +
    geom_hex(bins=70) +
    theme_bw() +
    scale_fill_continuous(type = "viridis", trans = "log10") +
    labs(title='Autosome variants, CADD>20', x='gnomAD fraction of alleles', y='EGA fraction of alleles') +
    stat_smooth(method='lm', color='red') +
    stat_cor(label.x.npc = 0.05, label.y.npc = 0.95)


# variants with the largest difference between classifications
vcf.20$gnomad.diff <- vcf.20$frac.alleles -vcf.20$gnomad.to.plot
vcf.20$gnomad.l2fc <- log2(vcf.20$frac.alleles / (vcf.20$gnomad.to.plot + 0.00001))
vcf.20$gene_change <- paste(vcf.20$hugoSymbol, vcf.20$proteinChange, sep='_')
vcf.20.meta <- vcf.20[,colnames(vcf.20)[!(colnames(vcf.20) %in% vcf.sample.names)]]
remove.cols <- c('secondaryVariantClassification',
                 'ncbiBuild',
                 'cDnaChange',
                 'codonChange')
vcf.20.meta <- vcf.20.meta[,!(colnames(vcf.20.meta) %in% remove.cols)]
# save these result tables
write.table(vcf.20, '~/pcloud_sync/project/EGA_data/results/germline/vcf_20.tsv', sep='\t', quote=F, row.names = T, col.names = T)
write.table(vcf.20.meta, '~/pcloud_sync/project/EGA_data/results/germline/vcf_20_meta.tsv', sep='\t', quote=F, row.names = T, col.names = T)

# top 100 different variants, high in EGA but low in gnomAD
vcf.20.topdiff <- vcf.20[vcf.20$gnomad.diff >= 0.05, ]
vcf.20.topdiff <- vcf.20.topdiff[order(vcf.20.topdiff$gnomad.l2fc, decreasing = T)[1:1000], colnames(vcf.20.topdiff)[!(colnames(vcf.20.topdiff) %in% vcf.sample.names)]]

remove.genes <- c('HYDIN')
vcf.20.topdiff <- vcf.20.topdiff[,!(colnames(vcf.20.topdiff) %in% remove.cols)]
vcf.20.topdiff <- vcf.20.topdiff[!(vcf.20.topdiff$hugoSymbol %in% remove.genes),]
# eliminate things found in >90% of patients because they're surely false positives
vcf.20.topdiff <- vcf.20.topdiff[vcf.20.topdiff$frac.pt.with.variant < 0.75, ]
# less thatn 20% in gnomad
vcf.20.topdiff <- vcf.20.topdiff[vcf.20.topdiff$gnomad.to.plot < 0.35, ]

# scatter of the top interesting variants
dim(vcf.20.topdiff)
ggplot(vcf.20.topdiff, aes(x=gnomad.to.plot, y=frac.alleles, label=hugoSymbol)) +
    geom_point(alpha=0.5) +
    theme_bw() +
    ylim(c(0,0.4)) +
    xlim(c(0,0.4)) +
    # geom_vline(xintercept = 0.36, col='firebrick') +
    # geom_hline(yintercept = 0.75, col='firebrick') +
    geom_abline(slope=1, col='firebrick') +
    labs(title='Autosome variants, CADD>20', x='gnomAD fraction of alleles', y='EGA fraction of alleles') +
    stat_dens2d_labels(geom = "text_repel", keep.fraction = 0.1)

#Add gene descriptions to table
library(biomaRt)
httr::set_config(httr::config(ssl_verifypeer = FALSE))
# change rownames of count table first
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene.ids <- unique(vcf.20.topdiff$hugoSymbol)
foo.desc <- getBM(attributes=c('external_gene_name',
                               'description'),
                  filters = 'external_gene_name',
                  values = gene.ids,
                  mart = ensembl)

gene.to.desc <- setNames(foo.desc$description, foo.desc$external_gene_name)
vcf.20.topdiff$gene_description <- sapply(gene.to.desc[vcf.20.topdiff$hugoSymbol], function(x) strsplit(x, split=" \\[")[[1]][1])

# some of these I've checked and they're actually common and uninteresting
common.rs.filter <- c('rs17667531', 'rs77375493')
vcf.20.topdiff <- vcf.20.topdiff[!(vcf.20.topdiff$dbSNP_ID %in% common.rs.filter),]

# save the table to share
write.table(vcf.20.topdiff, "~/pcloud_sync/project/EGA_data/results/germline/EGA_germline_mutations_top_gnomAD_difference.tsv", sep="\t", quote=F, row.names = F, col.names = T)


# some interesting rows, how many are homo/het
interesting.rows <- c('chr3_133356790_G_A_TOPBP1_p.S817L_22.7', 'chr16_88964557_C_G_CBFA2T3_p.R103P_21.6')
for (r in interesting.rows){
    vd <- vcf.20.data[r, ]
    print(r)
    print(table(as.character(vd[1,])))
}

# look for specific genes in the whole vcf table
vcf.bod <- vcf.20.meta[vcf.20.meta$hugoSymbol == "BOD1L1", ]
vcf.top <- vcf.20.meta[vcf.20.meta$hugoSymbol == "TOPBP1", ]
vcf.cbfa <- vcf.20.meta[vcf.20.meta$hugoSymbol == "CBFA2T3", ]
