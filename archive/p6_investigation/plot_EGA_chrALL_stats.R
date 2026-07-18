# repeating the plots about ROH in the EGA samples for all chromosomes
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)

# get the EGA sample annotations from the metadata file
ega_metadata <- read.table('~/local_data/EGA_data/ega_metadata_join_roh.tsv', sep='\t', header = T, comment.char = '', quote='')
rownames(ega_metadata) <- ega_metadata$sample
ega_loh <- read.table('~/local_data/EGA_data/ega_loh_metadata.tsv', sep='\t', header = T, comment.char = '', quote='')
rownames(ega_loh) <- ega_loh$patient

df_list <- lapply(1:22, function(i) {
    df <- read.table(paste0("~/local_data/EGA_data/roh_stats_df/chr", i,"_roh_stats.txt"), sep='\t', header=T, quote='')
    df <- df[,c(1,3,4)]
    colnames(df) <- c('sample', paste0('chr', i, '_roh_number'), paste0('chr', i, '_roh_total_length'))
    return(df)
})
df_list2 <- lapply(1:22, function(i) {
    df <- read.table(paste0("~/local_data/EGA_data/roh_stats_df/chr", i,"_roh_stats.txt"), sep='\t', header=T, quote='')
    df <- df[,c(1,3,4)]
    df$chr <- i
    return(df)
})

df_long <- do.call(rbind, df_list2)
df_long <- df_long[df_long$sample != 'PD6567b',]
# annotate with some info from metadata
df_long$tumor_normal <- ega_metadata[df_long$sample, "tumor_normal"]
df_long$study <- ega_metadata[df_long$sample, "study"]
df_long$patient <- ega_metadata[df_long$sample, "patient"]
df_long$chr9p_loh <- ega_loh[df_long$patient, "chr9p_loh"]
df_long$any_loh <- ega_loh[df_long$patient, "any_loh"]
df_long$any_chr_event <- ega_loh[df_long$patient, "any_chr_event"]

df_sum <- df_long[df_long$chr==1, ]
df_sum$roh_number <-aggregate(df_long$roh_number, by=list(df_long$sample), FUN=sum)$x
df_sum$roh_total_length <-aggregate(df_long$roh_total_length, by=list(df_long$sample), FUN=sum)$x

ggplot(df_long, aes(x=roh_number, y=roh_total_length)) +
    geom_point(alpha=0.25) + 
    theme_bw() + 
    labs(title='ROH: EGA data, all chromosomes, all samples')

ggplot(df_long, aes(x=roh_number, y=roh_total_length, color=study)) +
    geom_point(alpha=0.5) + 
    theme_bw() + 
    facet_wrap(.~chr)+ 
    labs(title='ROH: EGA data, per chromosome, all samples')

ggplot(df_long, aes(x=study, y=roh_total_length, fill=study)) +
    geom_boxplot() + 
    theme_bw() + 
    facet_wrap(.~chr, ncol=6)+ 
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    ylim(c(0,2.75e8)) + 
    labs(title='ROH: EGA data, per chromosome, all samples')

p1 <- ggplot(df_sum, aes(x=study, y=roh_total_length, fill=study)) +
    geom_violin(draw_quantiles = c(0.5)) + 
    theme_bw() + 
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    labs(title='ROH: EGA data, SUMMED CHROMOSOMES, all samples')
p2 <- ggplot(df_sum, aes(x=study, y=roh_number, fill=study)) +
    geom_violin(draw_quantiles = c(0.5)) + 
    theme_bw() + 
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    labs(title='ROH: EGA data, SUMMED CHROMOSOMES, all samples')
plot_grid(p1,p2, nrow=2)

# bimodal chr9 in PV
ggplot(df_long[df_long$chr==9,], aes(x=roh_total_length, fill=study)) +
    geom_histogram() + 
    theme_bw() + 
    facet_wrap(.~study)+ 
    labs(title='ROH: EGA data, chr9')
ggplot(df_long[df_long$chr==9,], aes(x=roh_number, fill=study)) +
    geom_histogram() + 
    theme_bw() + 
    facet_wrap(.~study)+ 
    labs(title='ROH: EGA data, chr9')

# startify by loh status
ggplot(df_long[!is.na(df_long$chr9p_loh) & 
               df_long$tumor_normal %in% c('N'),],
       aes(x=interaction(chr9p_loh, study), y=roh_total_length, fill=chr9p_loh)) +
    geom_boxplot() + 
    theme_bw() + 
    stat_compare_means(comparisons=list(c('FALSE.ET', 'TRUE.ET'),
                            c('FALSE.PMF', 'TRUE.PMF'),
                            c('FALSE.PV', 'TRUE.PV'))) + 
    facet_wrap(.~chr, ncol=6)+ 
    ylim(c(0,2.75e8)) +
    labs(title='ROH: EGA data, per chromosome, germline samples')+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggplot(df_long[!is.na(df_long$any_chr_event) & 
               df_long$tumor_normal %in% c('N'),],
       aes(x=interaction(any_chr_event, study), y=roh_total_length, fill=any_chr_event)) +
    geom_boxplot() + 
    theme_bw() + 
    stat_compare_means(comparisons=list(c('FALSE.ET', 'TRUE.ET'),
                            c('FALSE.PMF', 'TRUE.PMF'),
                            c('FALSE.PV', 'TRUE.PV'))) + 
    facet_wrap(.~chr, ncol=6)+ 
    ylim(c(0,2.75e8)) +
    labs(title='ROH: EGA data, per chromosome, germline samples') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p1 <- ggplot(df_sum[!is.na(df_sum$chr9p_loh) & 
               df_sum$tumor_normal %in% c('N'),],
       aes(x=interaction(chr9p_loh, study), y=roh_total_length, fill=chr9p_loh)) +
    geom_boxplot() + 
    theme_bw() + 
    stat_compare_means(comparisons=list(c('FALSE.ET', 'TRUE.ET'),
                            c('FALSE.PMF', 'TRUE.PMF'),
                            c('FALSE.PV', 'TRUE.PV'))) + 
    labs(title='ROH: EGA data, SUMMED CHROMOSOMES, germline samples') 
p2 <- ggplot(df_sum[!is.na(df_sum$any_chr_event) & 
               df_sum$tumor_normal %in% c('N'),],
       aes(x=interaction(any_chr_event, study), y=roh_total_length, fill=any_chr_event)) +
    geom_boxplot() + 
    theme_bw() + 
    stat_compare_means(comparisons=list(c('FALSE.ET', 'TRUE.ET'),
                            c('FALSE.PMF', 'TRUE.PMF'),
                            c('FALSE.PV', 'TRUE.PV'))) + 
    labs(title='ROH: EGA data, SUMMED CHROMOSOMES, germline samples') 
plot_grid(p1,p2, nrow=2)


df_merge <- data.frame(sample=df_list[[1]]$sample)
dfm2 <- do.call(cbind, lapply(df_list, function(x) x[,2:3]))
df_merge <- cbind(df_merge, dfm2)
# this sample has no data
df_merge[which(df_merge$chr1_roh_number==0),]
df_merge <- df_merge[df_merge$sample != 'PD6567b',]
df_melt <- melt(df_merge)

ggplot(df_merge, aes(x=chr1_roh_number, y=chr1_roh_total_length)) +
    geom_point()

# PD8943a really stands out as an outlier here
df_merge[df_merge$chr1_roh_number<100 & df_merge$chr1_roh_total_length>1.4e+08,]
