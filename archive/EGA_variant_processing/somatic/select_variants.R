# Consoldation of the EGA variants, further filtering for likely true somatic mutations.
# read necessary metadata
library(ggplot2)
library(reshape2)

metadata.sample.f <- '~/pcloud_sync/project/EGA_data/metadata/bam_metadata_uniq_samp.tsv'
metadata.patient.f <- '~/pcloud_sync/project/EGA_data/metadata/patient_metadata.tsv'
vcf.dir <- '~/local_data/EGA_data/02_variants_reference_exome/REDO_new_scores/filtered/'

metadata.sample <- read.table(metadata.sample.f, sep='\t', quote='', header=T, comment.char = '')
metadata.patient <- read.table(metadata.patient.f, sep='\t', quote='', header=T)
rownames(metadata.sample) <- metadata.sample$sample
rownames(metadata.patient) <- metadata.patient$patient

# OMIT A FEW SAMPLES WHICH ARE POORLY BEHAVED AND I BELIVE THEY'RE NOT PAIRS
# and something is weird and/or missing with patient PD8640
remove.patients <- c('PD6637', 'PD5003', 'PD6628', 'PD6652', 'PD8636', 'PD8640')
metadata.patient <- metadata.patient[!(metadata.patient$patient %in% remove.patients), ]
metadata.sample <- metadata.sample[!(metadata.sample$patient %in% remove.patients), ]

# only going to work with patients where we have solid tumor normal pairs
metadata.patient <- metadata.patient[metadata.patient$tumor.and.normal,]
metadata.sample <- metadata.sample[metadata.sample$patient %in% metadata.patient$patient, ]
table(metadata.patient$ega.class)
type.df <- data.frame(table(metadata.patient$ega.class))
type.df$frac <- type.df$Freq/sum(type.df$Freq)
type.df$x <- 'EGA'
type.df$Var1 <- factor(c('ET', 'PMF', 'PV'), levels=rev(c('ET', 'PV', 'PMF')))
ggplot(type.df, aes(fill=Var1, y=frac, x=x)) +
    geom_bar(position="stack", stat="identity") +
    coord_flip() +
    theme_bw() +
    scale_fill_brewer(palette = 'Set2', name='MPN') +
    labs(y='Fraction of patients', x='')

# for each patient, read in the VCF
vcf.list <- lapply(metadata.patient$patient, function(x) {
    print(x)
    read.table(file.path(vcf.dir, paste0(x, '.vcf')), sep='\t', quote='', header=T)
})
names(vcf.list) <- metadata.patient$patient

# implement some other quality filters
# MOST of these should already be applied
filter.min.T.AF <- 0.05
filter.min.CADD <- 20
filter.min.DP <- 10
filter.l2fc <- 1
filter.GATK.list <- c('PASS')
pt <- metadata.patient$patient[1]
vcf.list.filtered <- vcf.list
filter_vcf_df <- function(pt){
    samples <- sort(metadata.sample[metadata.sample$patient==pt, "sample"])
    sample.colnames <- c(sapply(samples, function(x) paste0(c('AF_', 'DP_'), x)))
    tumor.samples <- samples[metadata.sample[samples, "tumor_normal"] =='T']
    normal.samples <- samples[metadata.sample[samples, "tumor_normal"] =='N']
    tumor.sample.colnames <- c(sapply(tumor.samples, function(x) paste0(c('AF_', 'DP_'), x)))
    tumor.sample.af.colnames <- c(sapply(tumor.samples, function(x) paste0(c('AF_'), x)))
    tumor.sample.dp.colnames <- c(sapply(tumor.samples, function(x) paste0(c('DP_'), x)))
    normal.sample.colnames <- c(sapply(normal.samples, function(x) paste0(c('AF_', 'DP_'), x)))
    normal.sample.af.colnames <- c(sapply(normal.samples, function(x) paste0(c('AF_'), x)))
    normal.sample.dp.colnames <- c(sapply(normal.samples, function(x) paste0(c('DP_'), x)))
    vcf.df <- vcf.list[[pt]]

    # Tumor AF
    pass.min.T.AF <- apply(vcf.df[,tumor.sample.af.colnames, drop=F], 1, function(x) any(x >= filter.min.T.AF, na.rm = T))
    # Min CADD
    pass.min.CADD <- vcf.df$CADD_phred >= filter.min.CADD
    # min depth in at least one tumor and one normal
    pass.min.DP <- apply(vcf.df[,tumor.sample.dp.colnames, drop=F], 1, function(x) any(x >= filter.min.DP, na.rm = T)) &
                        apply(vcf.df[,normal.sample.dp.colnames, drop=F], 1, function(x) any(x >= filter.min.DP, na.rm = T))
    # filter on these three
    vcf.df <- vcf.df[pass.min.T.AF & pass.min.CADD & pass.min.DP, ]

    # min fold change
    # calculate that because its not a field yet
    e <- 0.001
    max.l2fc <- 2.5
    vcf.df$l2fc <- log2(vcf.df$mean.T.nna/(vcf.df$mean.N.nna + e))
    # cap at 2.5
    vcf.df$l2fc[vcf.df$l2fc > max.l2fc] <- max.l2fc
    # filter at threshold
    vcf.df <- vcf.df[vcf.df$l2fc >= filter.l2fc, ]

    # PASS GATK FILTERS
    # vcf.df <- vcf.df[vcf.df$FILTER %in% filter.GATK.list, ]

    return(vcf.df)
}

vcf.list.filtered <- lapply(metadata.patient$patient, filter_vcf_df)
names(vcf.list.filtered) <- metadata.patient$patient

# how many variants do we get per patient?
n.variants <- sapply(vcf.list.filtered, nrow)
names(n.variants) <- names(vcf.list.filtered)
hist(n.variants, breaks=20, main='Number of variants per patient', xlab='Number of variants', ylab='Number of patients')
head(sort(n.variants, decreasing = T))
# need to get to a format of one variant per line, indexed by patient and key

variant.distribution <- as.data.frame(sort(table(unlist(lapply(vcf.list.filtered, function(x) x$variantClassification))), decreasing = T))
variant.distribution$Var1 <- factor(variant.distribution$Var1, levels=variant.distribution$Var1)
ggplot(variant.distribution, aes(x=Var1, y=Freq)) +
    geom_bar(stat='identity')  +
    coord_flip() +
    theme_bw() +
    labs(x='Frequency', y='Variant type')

# combine this filtered list down to a single row per variant
# need to deal with diferent number of normals
# maybe take sum of depth and scaled average of af?
simplify_vcf_df <- function(pt){
    print(pt)
    vcf.df <- vcf.list.filtered[[pt]]
    vcf.df$patient <- pt
    vcf.df$MPN <- metadata.patient[pt, "ega.class"]
    samples <- sort(metadata.sample[metadata.sample$patient==pt, "sample"])
    sample.colnames <- c(sapply(samples, function(x) paste0(c('AF_', 'DP_'), x)))
    tumor.samples <- samples[metadata.sample[samples, "tumor_normal"] =='T']
    normal.samples <- samples[metadata.sample[samples, "tumor_normal"] =='N']
    tumor.sample.colnames <- c(sapply(tumor.samples, function(x) paste0(c('AF_', 'DP_'), x)))
    tumor.sample.af.colnames <- c(sapply(tumor.samples, function(x) paste0(c('AF_'), x)))
    tumor.sample.dp.colnames <- c(sapply(tumor.samples, function(x) paste0(c('DP_'), x)))
    normal.sample.colnames <- c(sapply(normal.samples, function(x) paste0(c('AF_', 'DP_'), x)))
    normal.sample.af.colnames <- c(sapply(normal.samples, function(x) paste0(c('AF_'), x)))
    normal.sample.dp.colnames <- c(sapply(normal.samples, function(x) paste0(c('DP_'), x)))

    keep.colnames <- c("patient","MPN","CHROM","POS","REF","ALT","FILTER","hugoSymbol","ncbiBuild","variantClassification","secondaryVariantClassification","variantType","proteinChange","genomeChange","cDnaChange","codonChange","CADD_phred","gene_protein","l2fc","score.new")
    vcf.df.simple <- vcf.df[,keep.colnames]
    vcf.df.simple$T.AF <- apply(vcf.df[,tumor.sample.af.colnames, drop=F], 1, function(x) mean(x, na.rm=T))
    vcf.df.simple$T.DP <- apply(vcf.df[,tumor.sample.dp.colnames, drop=F], 1, function(x) sum(x, na.rm=T))
    vcf.df.simple$N.AF <- apply(vcf.df[,normal.sample.af.colnames, drop=F], 1, function(x) mean(x, na.rm=T))
    vcf.df.simple$N.DP <- apply(vcf.df[,normal.sample.dp.colnames, drop=F], 1, function(x) sum(x, na.rm=T))
    return(vcf.df.simple)
}
vcf.list.filtered.simple <- lapply(metadata.patient$patient, simplify_vcf_df)
names(vcf.list.filtered.simple) <- metadata.patient$patient

vcf.simple.df <- do.call(rbind,vcf.list.filtered.simple)
head(vcf.simple.df)
vcf.simple.df <- vcf.simple.df[!is.na(vcf.simple.df$hugoSymbol),]

cs.vars <- vcf.simple.df[vcf.simple.df$hugoSymbol=="CS",]
calr.vars <- vcf.simple.df[vcf.simple.df$hugoSymbol=="CALR",]

# answer some basic questions about this
# what genes are most commonly mutated, and per MPN?
# subset to SNPs
keep.variant.types <- c('MISSENSE','NONSENSE')
v <- vcf.simple.df[vcf.simple.df$variantClassification %in% keep.variant.types,]
v$pt_gene <- paste(v$patient, v$hugoSymbol, sep='_')
v.uniq.pt.gene <- v[!(duplicated(v$pt_gene)),]
uniq.genes.count <- sort(table(v.uniq.pt.gene$hugoSymbol), decreasing = T)
uniq.genes.count.mpn <- table(v.uniq.pt.gene$hugoSymbol, v.uniq.pt.gene$MPN)
uniq.genes.count.mpn <- uniq.genes.count.mpn[order(rowSums(uniq.genes.count.mpn), decreasing = T),]
uniq.genes.freq.mpn <- t(apply(uniq.genes.count.mpn, 1, function(x) x / as.numeric(table(metadata.patient$ega.class))))


top.genes.mpn <- rownames(uniq.genes.count.mpn)[1:50]
ugc.df <- data.frame(uniq.genes.count)
ugf.df.mpn <- melt(uniq.genes.freq.mpn)

ggplot(ugc.df[1:50,], aes(x=Var1, y=Freq)) +
    geom_bar(stat='identity')  +
    coord_flip() +
    theme_bw() +
    labs(y='Number of patients with mutation', x='Gene',
         title='Commonly mutated genes: MISSENSE and NONSENSE variants')

ggplot(ugf.df.mpn[ugf.df.mpn$Var1%in% top.genes.mpn, ], aes(x=Var1, y=value)) +
    geom_bar(stat='identity')  +
    coord_flip() +
    theme_bw() +
    labs(y='Fraction of patients with mutation', x='Gene',
         title='Commonly mutated genes: MISSENSE and NONSENSE variants') +
    facet_wrap(.~Var2)


# write out this table with the filtered results
outf <- '~/pcloud_sync/project/EGA_data/results/combined_somatics_exome.tsv'
write.table(vcf.simple.df, outf, sep='\t', quote=F, row.names = F, col.names = T)
vcf.simple.df <- read.table('~/pcloud_sync/project/EGA_data/results/combined_somatics_exome.tsv', sep='\t', quote='', header=T)
