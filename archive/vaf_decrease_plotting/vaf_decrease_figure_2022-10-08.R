# load the VAF data from user 
# Sourced from this file on sharepoint: All_QC_CTP_201_20220914_updated_20221006.xlsx Tbl3 Fu VAF

library(ggplot2)
library(ggpubr)
library(cowplot)

df <- read.table('~/local_data/VAF_figure/All_QC_CTP_201_20220914_updated_20221006.tsv', sep='\t', quote='', header = T, fill=T, comment.char = '')
df$pt_gene <- paste(df$patient, df$Gene, sep='_')
sort(table(df[!duplicated(df$pt_gene), 'Gene']))
# threshold log10
max_abs_val <- 1
df$Log10Change_trunc <- df$Log10Change
df$Log10Change_trunc[df$Log10Change_trunc > max_abs_val] <- max_abs_val
df$Log10Change_trunc[df$Log10Change_trunc < -max_abs_val] <- -max_abs_val

# bucket the DaysPostDay0 into 5 buckets
df$day_bucket <- floor(df$DaysPostDay0 / 100)
df$Log2Change <- log2(10^df$Log10Change)
# it only makes sense to do this on the CALR and JAK2 genes - 12 and 17 patients with these data respectively. 
df_JAK2 <- df[df$Gene=='JAK2', ]
df_CALR <- df[df$Gene=='CALR', ]
df_MPL <- df[df$Gene=='MPL', ]

# a few manual changes
# ignore EOS when we have an EOT sample
df_JAK2 <- df_JAK2[!((df_JAK2$patient %in% c('057-201', '009-201', '021-201')) & (df_JAK2$timepoint=='EOS')), ]
df_CALR <- df_CALR[!((df_CALR$patient %in% c('050-202')) & (df_CALR$timepoint=='EOS')), ]
# remove 4th CALR sample of patient 056-201
df_CALR <- df_CALR[!((df_CALR$patient %in% c('056-201')) & (df_CALR$DaysPostDay0==420)), ]
# only one of the MPL muts
df_MPL <- df_MPL[df_MPL$key=='1_43815008_T_G', ]

df_JAK2 <- df_JAK2[order(df_JAK2$patient, df_JAK2$DaysPostDay0),]
df_JAK2$pt_sample_n <- do.call(c, sapply(table(df_JAK2$patient), function(x) 1:x))
df_JAK2$first_sample <- FALSE
df_JAK2$first_sample[which(!(duplicated(df_JAK2$patient)))] <- TRUE
df_JAK2$final_sample <- FALSE
df_JAK2$final_sample[which(!(duplicated(df_JAK2$patient,fromLast = T)))] <- TRUE
df_JAK2_final <- df_JAK2[df_JAK2$final_sample,]
df_JAK2_first <- df_JAK2[df_JAK2$first_sample,]
df_JAK2$pt_sample_n <- factor(df_JAK2$pt_sample_n, levels=c(1,2,3))

df_CALR <- df_CALR[order(df_CALR$patient, df_CALR$DaysPostDay0),]
df_CALR$pt_sample_n <- do.call(c, sapply(table(df_CALR$patient), function(x) 1:x))
df_CALR$first_sample <- FALSE
df_CALR$first_sample[which(!(duplicated(df_CALR$patient)))] <- TRUE
df_CALR$final_sample <- FALSE
df_CALR$final_sample[which(!(duplicated(df_CALR$patient,fromLast = T)))] <- TRUE
df_CALR_final <- df_CALR[df_CALR$final_sample,]
df_CALR_first <- df_CALR[df_CALR$first_sample,]
df_CALR$pt_sample_n <- factor(df_CALR$pt_sample_n, levels=c(1,2,3,4))

df_MPL$pt_sample_n <- factor(1, levels=c(1))

# define ordering of patients
df_JAK2$patient <- factor(df_JAK2$patient,
                          levels=unique(df_JAK2_first[order(df_JAK2_first$Log10Change), "patient"]))
df_CALR$patient <- factor(df_CALR$patient,
                          levels=unique(df_CALR_first[order(df_CALR_first$Log10Change), "patient"]))

# fill of the background box
df_JAK2$box_fill <- 'none'
df_JAK2[df_JAK2$final_sample, "box_fill"][df_JAK2[df_JAK2$first_sample, "Log10Change"] > 0] <- 'Increase'
df_JAK2[df_JAK2$final_sample, "box_fill"][df_JAK2[df_JAK2$first_sample, "Log10Change"] < 0] <- 'Decrease'
df_CALR$box_fill <- 'none'
df_CALR[df_CALR$final_sample, "box_fill"][df_CALR[df_CALR$first_sample, "Log10Change"] > 0] <- 'Increase'
df_CALR[df_CALR$final_sample, "box_fill"][df_CALR[df_CALR$first_sample, "Log10Change"] < 0] <- 'Decrease'
df_MPL$box_fill <- 'Decrease'

box_pal <- c(Increase= 'firebrick', Decrease='steelblue')
box_alpha <- 0.35
box_alpha_scale <- c(Increase= box_alpha, Decrease=box_alpha, none=0)
box_alpha_scale2 <- c(Increase= box_alpha, Decrease=box_alpha)


# Keeping multiple samples per patient
# But ignoring “time” component of X axis
p_JAK2 <- ggplot(df_JAK2, aes(x=as.character(DaysPostDay0), y=Log10Change_trunc, label = DaysPostDay0)) + 
    geom_rect(aes(ymin=-1.2, ymax=1.2, xmin=0.25, xmax=as.numeric(pt_sample_n)+0.75, fill=box_fill, alpha=box_fill), inherit.aes = F, show.legend = F) +
    scale_fill_manual(values=box_pal) +
    scale_alpha_manual(values=box_alpha_scale) +
    geom_bar(stat='identity',) +
    labs(x = 'Days post baseline', y='Log10 change \n VAF from baseline', title='JAK2 V617F') +
    facet_wrap(.~patient, nrow=1, scales='free_x') +
    theme_bw() +
    geom_hline(size=0.25, yintercept = 0) + 
    coord_cartesian(ylim=c(-1, 1)) +
    scale_x_discrete() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p_JAK2

p_CALR <- ggplot(df_CALR, aes(x=as.character(DaysPostDay0), y=Log10Change_trunc, label = DaysPostDay0)) + 
    geom_rect(aes(ymin=-1.2, ymax=1.2, xmin=0.25, xmax=as.numeric(pt_sample_n)+0.75, fill=box_fill, alpha=box_fill), inherit.aes = F, show.legend = F) +
    scale_fill_manual(values=box_pal) +
    scale_alpha_manual(values=box_alpha_scale) +
    geom_bar(stat='identity') +
    labs(x = 'Days post baseline', y='Log10 change \n VAF from baseline', title='CALR K385NCX or 34-52bp deletion') +
    facet_wrap(.~patient, nrow=1, scales = "free_x") +
    theme_bw() +
    geom_hline(size=0.25, yintercept = 0) + 
    coord_cartesian(ylim=c(-1, 1)) +
    scale_x_discrete() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p_CALR

p_MPL <- ggplot(df_MPL, aes(x=as.character(DaysPostDay0), y=Log10Change_trunc, label = DaysPostDay0)) + 
    geom_rect(aes(ymin=-1.2, ymax=1.2, xmin=0.25, xmax=as.numeric(pt_sample_n)+0.75, fill=box_fill, alpha=box_fill), inherit.aes = F, show.legend = T) +
    scale_fill_manual(values=box_pal) +
    scale_alpha_manual(values=box_alpha_scale2) +
    geom_bar(stat='identity') +
    labs(x = 'Days post baseline', y='Log10 change \n VAF from baseline', title='MPL W515G') +
    facet_wrap(.~patient, nrow=1, scales = "free_x") +
    theme_bw() +
    geom_hline(size=0.25, yintercept = 0) +
    coord_cartesian(ylim=c(-1, 1)) + 
    scale_x_discrete() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(fill='VAF change', alpha='VAF change')
p_MPL

p_combined <- plot_grid(p_JAK2, p_MPL, p_CALR, ncol=2, rel_widths = c(7,1))

# save them each as separate PDFs
pdf('~/local_data/VAF_figure/ASH_JAK2.pdf', height=2.75, width = 12); p_JAK2; dev.off()
png('~/local_data/VAF_figure/ASH_JAK2.png', height=2.75, width = 12, units='in', res = 300); p_JAK2; dev.off()
pdf('~/local_data/VAF_figure/ASH_CALR.pdf', height=2.75, width = 12); p_CALR; dev.off()
png('~/local_data/VAF_figure/ASH_CALR.png', height=2.75, width = 12, units='in', res = 300); p_CALR; dev.off()
pdf('~/local_data/VAF_figure/ASH_MPL.pdf', height=2.75, width = 3); p_MPL; dev.off()
png('~/local_data/VAF_figure/ASH_MPL.png', height=2.75, width = 3, units='in', res = 300); p_MPL; dev.off()

# make a legend to use in the plot


# plotting time on X axis - TRASH
ggplot(df_JAK2, aes(x=DaysPostDay0, y=Log10Change, label = DaysPostDay0)) + 
    geom_bar(stat='identity', width = 10) +
    labs(x = 'Patient followup sample', y='Log10 change in VAF from baseline', title='JAK2 V617F') +
    facet_wrap(.~patient, ncol=1,switch = 'y') +
    geom_hline(yintercept = 0, size=0.25) + 
    theme_bw() 


# looking at just a single timepoint in these data
tp <- 'ITPD169'
target_genes <- c('JAK2', 'CALR', 'MPL')
df_time <- df[df$timepoint==tp & df$Gene %in% target_genes, ]
df_time <- df_time[df_time$mut != 'MPL_W515S', ]
# order based on some gene, maybe JAK2 since it has the most information
jak2_order <- df_time[df_time$Gene=='JAK2', 'patient'][order(df_time[df_time$Gene=='JAK2', "Log10Change"])]
calr_order <- df_time[df_time$Gene=='CALR', 'patient'][order(df_time[df_time$Gene=='CALR', "Log10Change"])]
mpl_order <- df_time[df_time$Gene=='MPL', 'patient'][order(df_time[df_time$Gene=='MPL', "Log10Change"])]
df_time$patient <- factor(df_time$patient, levels = unique(c(jak2_order, calr_order, mpl_order)))

p1 <- ggplot(df_time, aes(x=Log10Change_trunc, y='')) + 
    geom_bar(stat='identity') +
    # labs(x = 'Patient followup sample', y='Log10 change \n VAF from baseline', title='JAK2 V617F') +
    facet_grid(rows=vars(patient), cols = vars(Gene)) +
    theme_bw() + 
    geom_vline(xintercept = 0, size=0.25)
    # geom_hline(size=0.25, yintercept = 0) + 
    # ylim(c(-1,1)) 
p1
