# association of germline variants with LOH, in the EGA data
library(ggplot2)
library(viridis)
library(gplots)
library(reshape2)

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


# copied from the LOH processing code
if(T){
metadata.sample.f <- '~/pcloud_sync/project/EGA_data/bam_metadata_uniq_samp.tsv'
metadata.patient.f <- '~/pcloud_sync/project/EGA_data/metadata/patient_metadata.tsv'
metadata.sample <- read.table(metadata.sample.f, sep='\t', quote='', header=T, comment.char = '')
metadata.patient <- read.table(metadata.patient.f, sep='\t', quote='', header=T)
rownames(metadata.sample) <- metadata.sample$sample
rownames(metadata.patient) <- metadata.patient$patient

# OMIT A FEW SAMPLES WHICH ARE POORLY BEHAVED AND I BELIVE THEY'RE NOT PAIRS
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

# association with specific germline variants
interesting.rows <- c('chr3_133356790_G_A_TOPBP1_p.S817L_22.7', 'chr16_88964557_C_G_CBFA2T3_p.R103P_21.6')
for (r in interesting.rows){
    metadata.patient[,r] <- as.integer(vcf.20[r, make.names(metadata.patient$patient)])
}

test.rows <- min(500, nrow(vcf.20.topdiff))
res.list <- lapply(rownames(vcf.20.topdiff)[1:test.rows], function(r) {
    print(r)
    m.test <- metadata.patient[,c("patient", "ega.class","any_event_count","loh_count","cnv_count", "chr9p_loh", 'chr8_amp',
                                  'chr1p_loh', 'any_event_but_chr9p')]
    m.test[,r] <- as.integer(vcf.20.topdiff[r, make.names(m.test$patient)])
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
rownames(res.df) <- vcf.20.topdiff$gene_change[1:test.rows]

res.df.plot <- -log10(res.df)
max.val <- 6
res.df.plot[res.df.plot >max.val] <- max.val
res.df.plot[res.df.plot <.001] <- 0
res.df.plot <- res.df.plot[order(apply(res.df.plot, 1, max), decreasing = T), ]

pdf('~/pcloud_sync/project/EGA_data/results/germline_loh_association_heatmap.pdf', height=10, width=10)
heatmap.2(res.df.plot[1:25, 1:28], Rowv = NA, Colv=NA, trace='none',
          col=viridis(32),
          margins = c(18,10),
          sepwidth = c(0.1,0.1),
          colsep = seq(7,ncol(res.df.plot), 7),
          key.title='-log10(p-value)',
          lhei = c(1,7))
dev.off()

top.hits <- rownames(res.df.plot)[1:10]
rdm <- melt(res.df.plot[,1:28], as.is = T)
rdm <- rdm[order(rdm$value ,decreasing = T),]
colnames(rdm) <- c('Gene_change', 'LOH', "-log10(p-value)")
rownames(rdm) <- NULL
rdm$p.adjust <- p.adjust(10^(-rdm$`-log10(p-value)`))
print(rdm[1:25,])

vcf.20.topdiff[vcf.20.topdiff$gene_change=='CNMD_p.V175I', ]
