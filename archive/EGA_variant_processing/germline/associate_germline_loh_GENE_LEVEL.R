# do the same association with the EGA data, germline variants and LOH,
# but this time I aggregate to the gene level first, rather than individual variants
# association of germline variants with LOH, in the EGA data
library(ggplot2)
library(cowplot)
library(viridis)
library(gplots)
library(reshape2)

vcf.20 <- read.table('~/pcloud_sync/project/EGA_data/results/germline/vcf_20.tsv', sep='\t', quote='', header = T)
remove.genes <- c('HYDIN')
# eliminate these certain genes
vcf.20 <- vcf.20[!(vcf.20$hugoSymbol %in% remove.genes),]
# eliminate variants in greater than 75% of patients
vcf.20 <- vcf.20[vcf.20$frac.pt.with.variant < 0.75, ]
# eliminate variants high in gnomAD which would cloud our aggregation with commmon crap
# this cutoff needs to be set fairly low
vcf.20 <- vcf.20[vcf.20$gnomad.to.plot < 0.1, ]

# copied from the LOH processing code
if(T){
metadata.sample.f <- '~/pcloud_sync/project/EGA_data/bam_metadata_uniq_samp.tsv'
metadata.patient.f <- '~/pcloud_sync/project/EGA_data/metadata/patient_metadata.tsv'
metadata.sample <- read.table(metadata.sample.f, sep='\t', quote='', header=T, comment.char = '')
metadata.patient <- read.table(metadata.patient.f, sep='\t', quote='', header=T)
rownames(metadata.sample) <- metadata.sample$sample
rownames(metadata.patient) <- metadata.patient$patient

# OMIT A FEW SAMPLES WHICH ARE POORLY BEHAVED AND I BELIVE THEY'RE NOT PAIRS
# and this one who doesnt have gemline calls yet?
remove.patients <- c('PD6637', 'PD5003', 'PD6628', 'PD6652', 'PD8636', 'PD8640', 'PD6647')
metadata.patient <- metadata.patient[!(metadata.patient$patient %in% remove.patients), ]
metadata.sample <- metadata.sample[!(metadata.sample$patient %in% remove.patients), ]

# only going to work with patients where we have solid tumor normal pairs
# metadata.patient <- metadata.patient[metadata.patient$tumor.and.normal,]
# metadata.sample <- metadata.sample[metadata.sample$patient %in% metadata.patient$patient, ]

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
}

# separate into matrices
vcf.20.data <- as.matrix(vcf.20[, make.names(metadata.patient$patient)])
vcf.20.meta <- vcf.20[, !(colnames(vcf.20) %in% make.names(metadata.patient$patient))]

# aggregation to the gene level..............
# for a patient, take the max of the call of the variants for that gene?
# maybe a big melt is the right way to do it
v2m <- melt(vcf.20.data, as.is = T)
v2m$gene <- rep(vcf.20.meta$hugoSymbol, times=ncol(vcf.20.data))

va <- aggregate.data.frame(v2m$value, by = list(v2m$Var2, v2m$gene), FUN=max)
colnames(va) <- c('patient', 'gene', 'genotype')
vam <- acast(va, gene~patient)

# number of patients with at least a het
vam.meta <- data.frame(gene = rownames(vam),
                       n.pt = apply(vam, 1, function(x) sum(x >0, na.rm = T)),
                       n.het = apply(vam, 1, function(x) sum(x ==1 , na.rm = T)),
                       n.homo = apply(vam, 1, function(x) sum(x ==2 , na.rm = T))
                           )
vam.meta$frac.pt <- vam.meta$n.pt / ncol(vam)
vam.meta$frac.het <- vam.meta$n.het / ncol(vam)
vam.meta$frac.homo <- vam.meta$n.homo / ncol(vam)

vam.meta <- vam.meta[order(vam.meta$frac.pt, decreasing = T),]

p1 <- ggplot(vam.meta, aes(x=frac.pt)) +
    geom_histogram() +
    theme_bw() +
    xlim(c(NA, 0.75)) +
    labs(x='Fraction of patients with variant in gene', y='Count of gene')
p2 <- ggplot(vam.meta, aes(x=frac.het)) +
    geom_histogram() +
    theme_bw() +
    xlim(c(NA, 0.75)) +
    labs(x='Fraction of patients with heterozygous variant in gene', y='Count of gene')
p3 <- ggplot(vam.meta, aes(x=frac.homo)) +
    geom_histogram() +
    theme_bw() +
    xlim(c(NA, 0.75)) +
    labs(x='Fraction of patients with homozygous variant in gene', y='Count of gene')

plot_grid(p1,p2,p3, nrow=1)


# do the big association with these genes
test.rows <- nrow(vam.meta)
res.list <- lapply(vam.meta$gene[1:test.rows], function(r) {
    print(r)
    m.test <- metadata.patient[,c("patient", "ega.class","any_event_count","loh_count","cnv_count", "chr9p_loh", 'chr8_amp',
                                  'chr1p_loh', 'any_event_but_chr9p')]
    m.test[,r] <- as.integer(vam[r, make.names(m.test$patient)])
    n.het <- sum(m.test[,r] == 1, na.rm = T)
    n.homo <- sum(m.test[,r] == 2, na.rm = T)
    # print(n.het)
    # print(n.homo)
    m.test$germline.any <- FALSE
    m.test$het <- FALSE
    m.test$homo <- FALSE
    m.test$germline.any[m.test[,r] > 0] <- TRUE
    m.test$het[m.test[,r] == 1] <- TRUE
    m.test$homo[m.test[,r] == 2] <- TRUE
    test.names <- c('ANY_EVENT', 'ANY_LOH', 'ANY_CNV', 'CHR9P_LOH', 'CHR8_AMP', 'CHR1P_LOH', 'ANY_EVENT_EXCEPT_9P')
    # do several fisher tests
    pvalues.all.germline.any <- tryCatch(c(fisher.test(table(m.test$germline.any, m.test$any_event_count>0))$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$loh_count>0))$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$cnv_count>0))$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$chr9p_loh))$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$chr8_amp))$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$chr1p_loh))$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$any_event_but_chr9p))$p.value
    ), error=function(e) return(rep(1,7)))
    names(pvalues.all.germline.any) <- paste0('germline_ANY_', test.names)
    pvalues.ET.germline.any <- tryCatch(c(fisher.test(table(m.test$germline.any, m.test$any_event_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$loh_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$cnv_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$chr9p_loh, m.test$ega.class)[,,"ET"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$chr8_amp, m.test$ega.class)[,,"ET"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$chr1p_loh, m.test$ega.class)[,,"ET"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$any_event_but_chr9p, m.test$ega.class)[,,"ET"])$p.value
    ), error=function(e) return(rep(1,7)))
    names(pvalues.ET.germline.any) <- paste0("germline_ANY_ET_", test.names)
    pvalues.PMF.germline.any <- tryCatch(c(fisher.test(table(m.test$germline.any, m.test$any_event_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$loh_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$cnv_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$chr9p_loh, m.test$ega.class)[,,"PMF"])$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$chr8_amp, m.test$ega.class)[,,"PMF"])$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$chr1p_loh, m.test$ega.class)[,,"PMF"])$p.value,
                                           fisher.test(table(m.test$germline.any, m.test$any_event_but_chr9p, m.test$ega.class)[,,"PMF"])$p.value
    ), error=function(e) return(rep(1,7)))
    names(pvalues.PMF.germline.any) <- paste0("germline_ANY_PMF_", test.names)
    pvalues.PV.germline.any <- tryCatch(c(fisher.test(table(m.test$germline.any, m.test$any_event_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$loh_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$cnv_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$chr9p_loh, m.test$ega.class)[,,"PV"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$chr8_amp, m.test$ega.class)[,,"PV"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$chr1p_loh, m.test$ega.class)[,,"PV"])$p.value,
                                          fisher.test(table(m.test$germline.any, m.test$any_event_but_chr9p, m.test$ega.class)[,,"PV"])$p.value
    ), error=function(e) return(rep(1,7)))
    names(pvalues.PV.germline.any) <- paste0("germline_ANY_PV_", test.names)

    # Heterozygotes
    pvalues.all.het <- tryCatch(c(fisher.test(table(m.test$het, m.test$any_event_count>0))$p.value,
                                  fisher.test(table(m.test$het, m.test$loh_count>0))$p.value,
                                  fisher.test(table(m.test$het, m.test$cnv_count>0))$p.value,
                                  fisher.test(table(m.test$het, m.test$chr9p_loh))$p.value,
                                  fisher.test(table(m.test$het, m.test$chr8_amp))$p.value,
                                  fisher.test(table(m.test$het, m.test$chr1p_loh))$p.value,
                                  fisher.test(table(m.test$het, m.test$any_event_but_chr9p))$p.value
    ), error=function(e) return(rep(1,7)))
    pvalues.ET.het <- tryCatch(c(fisher.test(table(m.test$het, m.test$any_event_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                 fisher.test(table(m.test$het, m.test$loh_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                 fisher.test(table(m.test$het, m.test$cnv_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                 fisher.test(table(m.test$het, m.test$chr9p_loh, m.test$ega.class)[,,"ET"])$p.value,
                                 fisher.test(table(m.test$het, m.test$chr8_amp, m.test$ega.class)[,,"ET"])$p.value,
                                 fisher.test(table(m.test$het, m.test$chr1p_loh, m.test$ega.class)[,,"ET"])$p.value,
                                 fisher.test(table(m.test$het, m.test$any_event_but_chr9p, m.test$ega.class)[,,"ET"])$p.value
    ), error=function(e) return(rep(1,7)))
    pvalues.PMF.het <- tryCatch(c(fisher.test(table(m.test$het, m.test$any_event_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                  fisher.test(table(m.test$het, m.test$loh_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                  fisher.test(table(m.test$het, m.test$cnv_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                  fisher.test(table(m.test$het, m.test$chr9p_loh, m.test$ega.class)[,,"PMF"])$p.value,
                                  fisher.test(table(m.test$het, m.test$chr8_amp, m.test$ega.class)[,,"PMF"])$p.value,
                                  fisher.test(table(m.test$het, m.test$chr1p_loh, m.test$ega.class)[,,"PMF"])$p.value,
                                  fisher.test(table(m.test$het, m.test$any_event_but_chr9p, m.test$ega.class)[,,"PMF"])$p.value
    ), error=function(e) return(rep(1,7)))
    pvalues.PV.het <- tryCatch(c(fisher.test(table(m.test$het, m.test$any_event_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                 fisher.test(table(m.test$het, m.test$loh_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                 fisher.test(table(m.test$het, m.test$cnv_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                 fisher.test(table(m.test$het, m.test$chr9p_loh, m.test$ega.class)[,,"PV"])$p.value,
                                 fisher.test(table(m.test$het, m.test$chr8_amp, m.test$ega.class)[,,"PV"])$p.value,
                                 fisher.test(table(m.test$het, m.test$chr1p_loh, m.test$ega.class)[,,"PV"])$p.value,
                                 fisher.test(table(m.test$het, m.test$any_event_but_chr9p, m.test$ega.class)[,,"PV"])$p.value
    ), error=function(e) return(rep(1,7)))
    names(pvalues.all.het) <- paste0('germline_HET_', test.names)
    names(pvalues.ET.het) <- paste0("germline_HET_ET_", test.names)
    names(pvalues.PMF.het) <- paste0("germline_HET_PMF_", test.names)
    names(pvalues.PV.het) <- paste0("germline_HET_PV_", test.names)

    # Homozygotes
    pvalues.all.homo <- tryCatch(c(fisher.test(table(m.test$homo, m.test$any_event_count>0))$p.value,
                                   fisher.test(table(m.test$homo, m.test$loh_count>0))$p.value,
                                   fisher.test(table(m.test$homo, m.test$cnv_count>0))$p.value,
                                   fisher.test(table(m.test$homo, m.test$chr9p_loh))$p.value,
                                   fisher.test(table(m.test$homo, m.test$chr8_amp))$p.value,
                                   fisher.test(table(m.test$homo, m.test$chr1p_loh))$p.value,
                                   fisher.test(table(m.test$homo, m.test$any_event_but_chr9p))$p.value
    ), error=function(e) return(rep(1,7)))
    pvalues.ET.homo <- tryCatch(c(fisher.test(table(m.test$homo, m.test$any_event_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$loh_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$cnv_count>0, m.test$ega.class)[,,"ET"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$chr9p_loh, m.test$ega.class)[,,"ET"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$chr8_amp, m.test$ega.class)[,,"ET"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$chr1p_loh, m.test$ega.class)[,,"ET"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$any_event_but_chr9p, m.test$ega.class)[,,"ET"])$p.value
    ), error=function(e) return(rep(1,7)))
    pvalues.PMF.homo <- tryCatch(c(fisher.test(table(m.test$homo, m.test$any_event_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                   fisher.test(table(m.test$homo, m.test$loh_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                   fisher.test(table(m.test$homo, m.test$cnv_count>0, m.test$ega.class)[,,"PMF"])$p.value,
                                   fisher.test(table(m.test$homo, m.test$chr9p_loh, m.test$ega.class)[,,"PMF"])$p.value,
                                   fisher.test(table(m.test$homo, m.test$chr8_amp, m.test$ega.class)[,,"PMF"])$p.value,
                                   fisher.test(table(m.test$homo, m.test$chr1p_loh, m.test$ega.class)[,,"PMF"])$p.value,
                                   fisher.test(table(m.test$homo, m.test$any_event_but_chr9p, m.test$ega.class)[,,"PMF"])$p.value
    ), error=function(e) return(rep(1,7)))
    pvalues.PV.homo <- tryCatch(c(fisher.test(table(m.test$homo, m.test$any_event_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$loh_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$cnv_count>0, m.test$ega.class)[,,"PV"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$chr9p_loh, m.test$ega.class)[,,"PV"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$chr8_amp, m.test$ega.class)[,,"PV"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$chr1p_loh, m.test$ega.class)[,,"PV"])$p.value,
                                  fisher.test(table(m.test$homo, m.test$any_event_but_chr9p, m.test$ega.class)[,,"PV"])$p.value
    ), error=function(e) return(rep(1,7)))

    names(pvalues.all.homo) <- paste0('germline_HOMO_', test.names)
    names(pvalues.ET.homo) <- paste0("germline_HOMO_ET_", test.names)
    names(pvalues.PMF.homo) <- paste0("germline_HOMO_PMF_", test.names)
    names(pvalues.PV.homo) <- paste0("germline_HOMO_PV_", test.names)
    return(c(pvalues.all.germline.any, pvalues.ET.germline.any, pvalues.PMF.germline.any, pvalues.PV.germline.any,
             pvalues.all.het, pvalues.ET.het, pvalues.PMF.het, pvalues.PV.het,
             pvalues.all.homo, pvalues.ET.homo, pvalues.PMF.homo, pvalues.PV.homo))
}
)


res.df <- do.call(rbind, res.list)
rownames(res.df) <- vam.meta$gene[1:test.rows]

res.df.plot <- -log10(res.df)
max.val <- 99
res.df.plot[res.df.plot >max.val] <- max.val
res.df.plot[res.df.plot <.001] <- 0
res.df.plot <- res.df.plot[order(apply(res.df.plot, 1, max), decreasing = T), ]
res.df.plot[1:5,1:5]

pdf('~/pcloud_sync/project/EGA_data/results/germline_loh_association_heatmap_GENE_LEVEL.pdf', height=10, width=10)
heatmap.2(res.df.plot[1:20, 1:28], Rowv = NA, Colv=NA, trace='none',
          col=viridis(32),
          margins = c(18,6),
          sepwidth = c(0.1,0.1),
          colsep = seq(7,ncol(res.df.plot), 7),
          key.title='-log10(p-value)',
          lhei = c(1,7))
dev.off()

top.hits <- rownames(res.df.plot)[1:10]
vam.meta[top.hits, ]
rdm <- melt(res.df.plot[,1:28], as.is = T)
rdm <- rdm[order(rdm$value ,decreasing = T),]
colnames(rdm) <- c('Gene', 'LOH', "-log10(p-value)")
rownames(rdm) <- NULL
print(rdm[1:25,])
rdm$p.adjust <- p.adjust(10^(-rdm$`-log10(p-value)`))
