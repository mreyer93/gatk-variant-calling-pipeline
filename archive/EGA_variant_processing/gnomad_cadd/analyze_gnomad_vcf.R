# Analysis of the CADD-scored gnomad data
# right now for BOD1L1 and MAP1B
# which have been subsetted to the AML targes for accurate comparison with our data
library(ggplot2)

############################################################################################################
## BOD1L1 ##################################################################################################
############################################################################################################
gene <- 'BOD1L1'
og.gnomad <- read.table('~/data/gnomad_CADD/gnomAD_BOD1L1.csv', sep=',', header=T)
og.gnomad$key <- paste(og.gnomad$Chromosome, og.gnomad$Position, og.gnomad$Reference, og.gnomad$Alternate, sep='_')
# output for tables and figures
outdir <- '~/pcloud_sync/project/gnomad_germline/BOD1L1'

# cadd scored version
cadd.df <- read.table('~/data/gnomad_CADD/gnomAD_BOD1L1_cadd_subset.vcf', sep='\t', quote='')
colnames(cadd.df) <- c("Chromosome", "Position", "Reference", "Alternate", "Type", "Length", "ConsDetail", "GeneName", "PHRED")
cadd.df$key <- paste(cadd.df$Chromosome, cadd.df$Position, cadd.df$Reference, cadd.df$Alternate, sep='_')
# verify keys match up
all(cadd.df$key %in% og.gnomad$key)

# there are multiple annotations per positions
# take only the first for these figures
cadd.df.uniq <- cadd.df[!duplicated(cadd.df$key),]
rownames(cadd.df.uniq) <- cadd.df.uniq$key

# histogram of CADD scores
pdf(file.path(outdir, 'gnomAD_CADD_histogram.pdf'))
ggplot(cadd.df.uniq, aes(x=PHRED)) +
    geom_histogram(binwidth = 1) +
    theme_bw() +
    geom_vline(xintercept = 20, col='firebrick', lty=2) +
    labs(title=paste('gnomAD variants:', gene, '(Subsetted to IDT AML Targets)'), subtitle = 'CADD-score distribution')
dev.off()

# merge into the gnomad table
og.gnomad.sub <- og.gnomad[og.gnomad$key %in% cadd.df.uniq$key, ]
og.gnomad.sub$Source <- cadd.df.uniq[og.gnomad.sub$key, "PHRED"]
colnames(og.gnomad.sub)[6] <- 'PHRED'

pdf(file.path(outdir, 'gnomAD_CADD_vs_position.pdf'))
ggplot(og.gnomad.sub[order(og.gnomad.sub$Allele.Frequency),], aes(x=Position, y=PHRED, col= log10(Allele.Frequency))) +
    geom_point() +
    scale_color_viridis_c() +
    theme_bw() +
    geom_hline(yintercept = 20, col='firebrick', lty=2) +
    labs(title=paste('gnomAD variants:', gene, '(Subsetted to IDT AML Targets)'), subtitle = 'CADD-score vs Position, colored by allele frequency')
dev.off()

pdf(file.path(outdir, 'gnomAD_CADD_vs_AF.pdf'))
ggplot(og.gnomad.sub, aes(x=log10(Allele.Frequency), y=PHRED)) +
    geom_point(alpha=0.5) +
    theme_bw() +
    geom_smooth(method='lm', formula= y~x, col='darkblue') +
    labs(title=paste('gnomAD variants:', gene, '(Subsetted to IDT AML Targets)'), subtitle = 'CADD-score vs Population allele frequency')
dev.off()

# linear fit of CADD and log allele frequency
og.gnomad.sub$af.log <- log10(og.gnomad.sub$Allele.Frequency)
fit <- lm(PHRED~af.log, data=og.gnomad.sub)
summary(fit)

# some population-AF scaled density metric..... how?
View(og.gnomad.sub[order(og.gnomad.sub$PHRED, decreasing = T), ])

gnomad.g20 <- og.gnomad.sub[og.gnomad.sub$PHRED >=20, ]
View(gnomad.g20[order(gnomad.g20$Allele.Frequency, decreasing = T), ])

sort(table(og.gnomad.sub$VEP.Annotation), decreasing = T)
sort(table(gnomad.g20$VEP.Annotation), decreasing = T)


pdf(file.path(outdir, 'gnomAD_CADD20_AF_histogram.pdf'))
ggplot(gnomad.g20, aes(x=af.log)) +
    geom_histogram(binwidth = 0.2) +
    theme_bw() +
    labs(title=paste('gnomAD variants:', gene, '(Subsetted to IDT AML Targets)'), subtitle = 'AF distribution of CADD>=20 variants',
         x='log allele frequency')
dev.off()

write.table(gnomad.g20, '~/pcloud_sync/project/gnomad_germline/BOD1L1/BOD1L1_gnomad_g20.tsv', sep = '\t', quote=F, row.names = F, col.names = T)
write.table(og.gnomad.sub, '~/pcloud_sync/project/gnomad_germline/BOD1L1/BOD1L1_gnomad_original_sub.tsv', sep = '\t', quote=F, row.names = F, col.names = T)

############################################################################################################
## MAP1B ###################################################################################################
############################################################################################################
gene <- 'MAP1B'
og.gnomad <- read.table('~/pcloud_sync/project/gnomad_germline/MAP1B/gnomAD_MAP1B.csv', sep=',', header=T)
og.gnomad$key <- paste(og.gnomad$Chromosome, og.gnomad$Position, og.gnomad$Reference, og.gnomad$Alternate, sep='_')
# output for tables and figures
outdir <- '~/pcloud_sync/project/gnomad_germline/MAP1B'

# do restriction to AML targets
aml.targets <- read.table('~/local_data/references/AML_261_IDT/targets_w100_hg19.st.bed', sep='\t', header=F)
aml.targets.chr5 <- aml.targets[aml.targets$V1 == 'chr5', ]

# oh look they're all in the targets :)
og.gnomad$in.target <- sapply(og.gnomad$Position, function(x) any(aml.targets$V2 <=x & aml.targets$V3 >=x) )
table(og.gnomad$in.target)

# cadd scored version
cadd.df <- read.table('~/pcloud_sync/project/gnomad_germline/MAP1B/gnomAD_MAP1B_cadd_scored_simple.tsv', sep='\t', quote='')
colnames(cadd.df) <- c("Chromosome", "Position", "Reference", "Alternate", "Type", "Length", "ConsDetail", "GeneName", "PHRED")
cadd.df$key <- paste(cadd.df$Chromosome, cadd.df$Position, cadd.df$Reference, cadd.df$Alternate, sep='_')
# verify keys match up
all(cadd.df$key %in% og.gnomad$key)

# there are multiple annotations per positions
# take only the first for these figures
cadd.df.uniq <- cadd.df[!duplicated(cadd.df$key),]
rownames(cadd.df.uniq) <- cadd.df.uniq$key

# histogram of CADD scores
pdf(file.path(outdir, 'gnomAD_CADD_histogram.pdf'))
ggplot(cadd.df.uniq, aes(x=PHRED)) +
    geom_histogram(binwidth = 1) +
    theme_bw() +
    geom_vline(xintercept = 20, col='firebrick', lty=2) +
    labs(title=paste('gnomAD variants:', gene, '(Subsetted to IDT AML Targets)'), subtitle = 'CADD-score distribution')
dev.off()

# merge into the gnomad table
og.gnomad.sub <- og.gnomad[og.gnomad$key %in% cadd.df.uniq$key, ]
og.gnomad.sub$Source <- cadd.df.uniq[og.gnomad.sub$key, "PHRED"]
colnames(og.gnomad.sub)[6] <- 'PHRED'

pdf(file.path(outdir, 'gnomAD_CADD_vs_position.pdf'))
ggplot(og.gnomad.sub[order(og.gnomad.sub$Allele.Frequency),], aes(x=Position, y=PHRED, col= log10(Allele.Frequency))) +
    geom_point() +
    scale_color_viridis_c() +
    theme_bw() +
    geom_hline(yintercept = 20, col='firebrick', lty=2) +
    labs(title=paste('gnomAD variants:', gene, '(Subsetted to IDT AML Targets)'), subtitle = 'CADD-score vs Position, colored by allele frequency')
dev.off()

pdf(file.path(outdir, 'gnomAD_CADD_vs_AF.pdf'))
ggplot(og.gnomad.sub, aes(x=log10(Allele.Frequency), y=PHRED)) +
    geom_point(alpha=0.5) +
    theme_bw() +
    geom_smooth(method='lm', formula= y~x, col='darkblue') +
    labs(title=paste('gnomAD variants:', gene, '(Subsetted to IDT AML Targets)'), subtitle = 'CADD-score vs Population allele frequency')
dev.off()

# linear fit of CADD and log allele frequency
og.gnomad.sub$af.log <- log10(og.gnomad.sub$Allele.Frequency)
fit <- lm(PHRED~af.log, data=og.gnomad.sub)
summary(fit)

# some population-AF scaled density metric..... how?
View(og.gnomad.sub[order(og.gnomad.sub$PHRED, decreasing = T), ])

gnomad.g20 <- og.gnomad.sub[og.gnomad.sub$PHRED >=20, ]
View(gnomad.g20[order(gnomad.g20$Allele.Frequency, decreasing = T), ])

sort(table(og.gnomad.sub$VEP.Annotation), decreasing = T)
sort(table(gnomad.g20$VEP.Annotation), decreasing = T)


pdf(file.path(outdir, 'gnomAD_CADD20_AF_histogram.pdf'))
ggplot(gnomad.g20, aes(x=af.log)) +
    geom_histogram(binwidth = 0.2) +
    theme_bw() +
    labs(title=paste('gnomAD variants:', gene, '(Subsetted to IDT AML Targets)'), subtitle = 'AF distribution of CADD>=20 variants',
         x='log allele frequency')
dev.off()

write.table(gnomad.g20, '~/pcloud_sync/project/gnomad_germline/MAP1B/MAP1B_gnomad_g20.tsv', sep = '\t', quote=F, row.names = F, col.names = T)
write.table(og.gnomad.sub, '~/pcloud_sync/project/gnomad_germline/MAP1B/MAP1B_gnomad_original_sub.tsv', sep = '\t', quote=F, row.names = F, col.names = T)

# what's the rate among the population?
# assuming one variant per individual?

total.variant.alleles <- sum(gnomad.g20$Allele.Count)
median.alleles <- median(gnomad.g20$Allele.Number)

total.variant.alleles / (median.alleles /2) * 100
