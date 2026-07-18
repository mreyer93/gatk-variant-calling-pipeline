# Consoldation of the EGA variants, further filtering for likely true somatic mutations.
# read necessary metadata
library(ggplot2)
library(viridis)
library(gplots)
library(reshape2)

metadata.sample.f <- '~/pcloud_sync/project/EGA_data/bam_metadata_uniq_samp.tsv'
metadata.patient.f <- '~/pcloud_sync/project/EGA_data/metadata/patient_metadata.tsv'
vcf.dir <- '~/local_data/EGA_data/02_variants_reference_exome/REDO_new_scores/filtered_GATK/'
vcf.simple.f <- '~/pcloud_sync/project/EGA_data/results/combined_somatics_exome.tsv'

metadata.sample <- read.table(metadata.sample.f, sep='\t', quote='', header=T, comment.char = '')
metadata.patient <- read.table(metadata.patient.f, sep='\t', quote='', header=T)
rownames(metadata.sample) <- metadata.sample$sample
rownames(metadata.patient) <- metadata.patient$patient
vcf.simple.df <-read.table(vcf.simple.f, sep='\t', quote='', header=T)

# OMIT A FEW SAMPLES WHICH ARE POORLY BEHAVED AND I BELIVE THEY'RE NOT PAIRS
remove.patients <- c('PD6637', 'PD5003', 'PD6628', 'PD6652', 'PD8636')
metadata.patient <- metadata.patient[!(metadata.patient$patient %in% remove.patients), ]
metadata.sample <- metadata.sample[!(metadata.sample$patient %in% remove.patients), ]

# only going to work with patients where we have solid tumor normal pairs
metadata.patient <- metadata.patient[metadata.patient$tumor.and.normal,]
metadata.sample <- metadata.sample[metadata.sample$patient %in% metadata.patient$patient, ]

# LOH events
loh.f <- '~/pcloud_sync/project/EGA_data/results/loh_events.tsv'
loh.df <- read.table(loh.f, sep='\t', quote='', header=T)
loh <- loh.df[, c(1, 9,10,11,12)]
# remove PD4060 as we don't have variant calls for this patient yet
loh <- loh[loh$patient != 'PD4060', ]

loh.m <- melt(loh, measure.vars = c('event.1', 'event.2', 'event.3'))
loh.m <- loh.m[loh.m$value != '', ]
colnames(loh.m)[4] <- 'event'
loh.m <- loh.m[, c("patient", "ega.class", "event")]
loh.m <- loh.m[order(loh.m$patient), ]
# no none here
loh.m.nn <- loh.m[loh.m$event !='none', ]

# remove PD4060 as we don't have variant calls for this patient yet
metadata.patient <- metadata.patient[metadata.patient$patient != 'PD4060', ]


# some qc on my calls
# make weak the same as normal
loh.m.nn$event.noweak <- gsub('weak ', '', loh.m.nn$event)
loh.m.nn$is.weak <- sapply(loh.m.nn$event, function(x) length(grep('weak', x))>0)
loh.m.nn$is.weak[loh.m.nn$is.weak==''] <- FALSE
loh.m.nn$chr_arm <- sapply(loh.m.nn$event.noweak, function(x) strsplit(x, split=' ')[[1]][1])
loh.m.nn$chr <- sapply(loh.m.nn$chr_arm, function(x) gsub('p', '', gsub('q', '', gsub('mid', '',x))))
loh.m.nn$event.type <- sapply(loh.m.nn$event.noweak, function(x) strsplit(x, split=' ')[[1]][2])

# table of LOH events by MPN
# any event
metadata.patient$any_event_count <- 0
metadata.patient[names(table(loh.m.nn$patient)), 'any_event_count'] <- table(loh.m.nn$patient)
metadata.patient$any_event_count_noweak <- 0
metadata.patient[names(table(loh.m.nn[!loh.m.nn$is.weak,"patient"])), 'any_event_count_noweak'] <- table(loh.m.nn[!loh.m.nn$is.weak,"patient"])
# specific LOH
metadata.patient$loh_count <- 0
metadata.patient[names(table(loh.m.nn[loh.m.nn$event.type=='loh',"patient"])), 'loh_count'] <- table(loh.m.nn[loh.m.nn$event.type=='loh',"patient"])
metadata.patient$loh_count_noweak <- 0
metadata.patient[names(table(loh.m.nn[loh.m.nn$event.type=='loh' & !loh.m.nn$is.weak ,"patient"])), 'loh_count_noweak'] <- table(loh.m.nn[loh.m.nn$event.type=='loh' &!loh.m.nn$is.weak,"patient"])
# specific to CNV
metadata.patient$cnv_count <- 0
metadata.patient[names(table(loh.m.nn[loh.m.nn$event.type %in% c('cn_amp', 'cn_loss'),"patient"])), 'cnv_count'] <- table(loh.m.nn[loh.m.nn$event.type %in% c('cn_amp', 'cn_loss'),"patient"])
metadata.patient$cnv_count_noweak <- 0
metadata.patient[names(table(loh.m.nn[loh.m.nn$event.type %in% c('cn_amp', 'cn_loss') & !loh.m.nn$is.weak ,"patient"])), 'cnv_count_noweak'] <- table(loh.m.nn[loh.m.nn$event.type %in% c('cn_amp', 'cn_loss') &!loh.m.nn$is.weak,"patient"])
# chr9 LOH
chr9.events <- c('chr9p loh', 'chr9p cn_amp', 'chr9 cn_amp')
metadata.patient$chr9p_loh <- FALSE
metadata.patient[loh.m.nn$patient[loh.m.nn$event.noweak %in%chr9.events], 'chr9p_loh'] <- TRUE
# chr8 cn_amp
metadata.patient$chr8_amp <- FALSE
metadata.patient[loh.m.nn$patient[loh.m.nn$event.noweak == 'chr8 cn_amp'], 'chr8_amp'] <- TRUE
#chr1p loh
metadata.patient$chr1p_loh <- FALSE
metadata.patient[loh.m.nn$patient[loh.m.nn$event.noweak == 'chr1p loh'], 'chr1p_loh'] <- TRUE
#ANYTHING excexpt chr9p
metadata.patient$any_event_but_chr9p <- FALSE
metadata.patient[loh.m.nn$patient[!(loh.m.nn$event.noweak %in%chr9.events)], 'any_event_but_chr9p'] <- TRUE

ggplot(metadata.patient, aes(x=loh_count)) +
    geom_bar(aes(y = ..prop..), stat="count") +
    facet_grid(.~ega.class) +
    theme_bw() +
    labs(x='Number of LOH or CNV events', y='Fraction of patients in group')

# Which events are the most common
loh.tab <- data.frame(table(loh.m.nn$event.noweak))
loh.tab <- loh.tab[order(loh.tab$Freq, decreasing = T), ]
loh.tab$Var1 <- factor(loh.tab$Var1, levels=loh.tab$Var1)
ggplot(loh.tab, aes(x=Var1, y=Freq))+
    geom_bar(stat='identity') +
    coord_flip() +
    theme_bw() +
    labs(x='Event', y='Count')

# make this plot by EGA class
ggplot(loh.m.nn, aes(x=event.noweak))+
    geom_bar(stat='count') +
    coord_flip() +
    theme_bw() +
    labs(x='Event', y='Count') +
    facet_wrap(.~ega.class)


# variant types to test
keep.variant.types <- c('MISSENSE','NONSENSE')
v <- vcf.simple.df[vcf.simple.df$variantClassification %in% keep.variant.types,]
v$pt_gene <- paste(v$patient, v$hugoSymbol, sep='_')
v.uniq.pt.gene <- v[!(duplicated(v$pt_gene)),]
uniq.genes.count <- sort(table(v.uniq.pt.gene$hugoSymbol), decreasing = T)
uniq.genes.count.mpn <- table(v.uniq.pt.gene$hugoSymbol, v.uniq.pt.gene$MPN)
uniq.genes.count.mpn <- uniq.genes.count.mpn[order(rowSums(uniq.genes.count.mpn), decreasing = T),]

test.genes <- names(uniq.genes.count)[1:1000]
test.gene.chr <- sapply(test.genes, function(g) vcf.simple.df[vcf.simple.df$hugoSymbol==g, "CHROM"][1])
res.list <- lapply(test.genes, function(g) {
    print(g)
    m.test <- metadata.patient[,c("patient", "ega.class","any_event_count","loh_count","cnv_count", "chr9p_loh", 'chr8_amp',
                                  'chr1p_loh', 'any_event_but_chr9p')]
    m.test$mutation <- FALSE
    m.test$mutation[m.test$patient %in% unique(vcf.simple.df[vcf.simple.df$hugoSymbol==g, "patient"])] <- TRUE
    test.names <- c('ANY_EVENT', 'ANY_LOH', 'ANY_CNV', 'CHR9P_LOH', 'CHR8_AMP', 'CHR1P_LOH', 'ANY_EVENT_EXCEPT_9P')
    # do several fisher tests
    pvalues.all <- c(fisher.test(table(m.test$mutation, m.test$any_event_count>0))$p.value,
                 fisher.test(table(m.test$mutation, m.test$loh_count>0))$p.value,
                 fisher.test(table(m.test$mutation, m.test$cnv_count>0))$p.value,
                 fisher.test(table(m.test$mutation, m.test$chr9p_loh))$p.value,
                 fisher.test(table(m.test$mutation, m.test$chr8_amp))$p.value,
                 fisher.test(table(m.test$mutation, m.test$chr1p_loh))$p.value,
                 fisher.test(table(m.test$mutation, m.test$any_event_but_chr9p))$p.value
    )
    names(pvalues.all) <- test.names
    pvalues.ET <- c(fisher.test(table(m.test$mutation, m.test$any_event_count>0, m.test$ega.class)[,,"ET"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$loh_count>0, m.test$ega.class)[,,"ET"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$cnv_count>0, m.test$ega.class)[,,"ET"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$chr9p_loh, m.test$ega.class)[,,"ET"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$chr8_amp, m.test$ega.class)[,,"ET"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$chr1p_loh, m.test$ega.class)[,,"ET"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$any_event_but_chr9p, m.test$ega.class)[,,"ET"])$p.value
    )
    names(pvalues.ET) <- paste0("ET_", test.names)
    pvalues.PMF <- c(fisher.test(table(m.test$mutation, m.test$any_event_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                     fisher.test(table(m.test$mutation, m.test$loh_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                     fisher.test(table(m.test$mutation, m.test$cnv_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                     fisher.test(table(m.test$mutation, m.test$chr9p_loh, m.test$ega.class)[,,"PMF"])$p.value,
                     fisher.test(table(m.test$mutation, m.test$chr8_amp, m.test$ega.class)[,,"PMF"])$p.value,
                     fisher.test(table(m.test$mutation, m.test$chr1p_loh, m.test$ega.class)[,,"PMF"])$p.value,
                     fisher.test(table(m.test$mutation, m.test$any_event_but_chr9p, m.test$ega.class)[,,"PMF"])$p.value
    )
    names(pvalues.PMF) <- paste0("PMF_", test.names)
    pvalues.PV <- c(fisher.test(table(m.test$mutation, m.test$any_event_count>0, m.test$ega.class)[,,"PV"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$loh_count>0, m.test$ega.class)[,,"PV"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$cnv_count>0, m.test$ega.class)[,,"PV"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$chr9p_loh, m.test$ega.class)[,,"PV"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$chr8_amp, m.test$ega.class)[,,"PV"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$chr1p_loh, m.test$ega.class)[,,"PV"])$p.value,
                    fisher.test(table(m.test$mutation, m.test$any_event_but_chr9p, m.test$ega.class)[,,"PV"])$p.value
    )
    names(pvalues.PV) <- paste0("PV_", test.names)
    return(c(pvalues.all, pvalues.ET, pvalues.PMF, pvalues.PV))
}
)

res.df <- do.call(rbind, res.list)
rownames(res.df) <- test.genes

res.df.plot <- -log10(res.df)
max.val <- 4
res.df.plot[res.df.plot >max.val] <- max.val
res.df.plot[res.df.plot <.001] <- 0
res.df.plot <- res.df.plot[order(apply(res.df.plot, 1, max), decreasing = T), ]

pdf('~/pcloud_sync/project/EGA_data/results/loh_mutation_association_heatmap.pdf', height=12, width=8)
heatmap.2(res.df.plot[1:50,], Rowv = NA, Colv=NA, trace='none',
          col=viridis(32),
          margins = c(13,6),
          colsep=c(7, 14,21,28),
          key.title='-log10(p-value)',
          lhei = c(1,7))
dev.off()

# aggregate to pick a list of top gene-disease pairs
rm <- melt(res.df)
rm <- rm[order(rm$value, decreasing = F),]
min.p <- 0.05
min.val <- -log10(min.p)
rm <- rm[rm$value <= min.p, ]
rm$value <- signif(rm$value, 2)
colnames(rm) <- c('gene', 'group', 'p')
rm$gene <- as.character(rm$gene)
rm$group <- as.character(rm$group)
rm$p <- as.numeric(rm$p)
rm$chr <- test.gene.chr[rm$gene]

# get total number of mutated patients
rm$mutated.patients <- uniq.genes.count[rm$gene]
rm$mutated.patients.ET <- uniq.genes.count.mpn[rm$gene, "ET"]
rm$mutated.patients.PMF <- uniq.genes.count.mpn[rm$gene, "PMF"]
rm$mutated.patients.PV <- uniq.genes.count.mpn[rm$gene, "PV"]

# write out
outf <- '~/pcloud_sync/project/EGA_data/results/loh_mutation_association.tsv'
write.table(rm, outf, sep='\t', quote=F, row.names = F, col.names = T)

test.pts <- c("PD5848","PD6568","PD7279","PD8625","PD8625")
loh.m.nn[loh.m.nn$patient %in% test.pts,]
