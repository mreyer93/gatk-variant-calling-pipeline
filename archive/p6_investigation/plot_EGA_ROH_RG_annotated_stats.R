# repeating the plots about ROH in the EGA samples for all chromosomes
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggExtra)

# get the EGA sample annotations from the metadata file
ega_metadata <- read.table('~/local_data/EGA_data/ega_metadata_join_roh.tsv', sep='\t', header = T, comment.char = '', quote='')
rownames(ega_metadata) <- ega_metadata$sample
ega_loh <- read.table('~/local_data/EGA_data/ega_loh_metadata.tsv', sep='\t', header = T, comment.char = '', quote='')
rownames(ega_loh) <- ega_loh$patient

df_list <- lapply(1:22, function(i) {
    df <- read.table(paste0("/path/to/data/p6/EGA_haplotypecaller/roh_annotated/roh_chr", i,"_RG.txt"), sep='\t', header=T, quote='')
    colnames(df) <-  c('class', 'sample', 'chr', 'start', 'end', 'length', 'num_markers', 'quality')
    return(df)
})
df <- do.call(rbind, df_list)
df <- df[df$sample!='PD6567b', ]
p <- ggplot(df, aes(x=quality, y=length)) + 
    geom_point(alpha=0.1) + 
    scale_y_log10()
p <- ggMarginal(p, type='histogram')
p
df[order(df$length, decreasing = T), ][1:15, ]

# if working with just chr9
# df <- read.table('/path/to/data/p6/EGA_haplotypecaller/roh_annotated/roh_chr9_RG.txt', sep='\t', quote='')
# colnames(df) <-  c('class', 'sample', 'chr', 'start', 'end', 'length', 'num_markers', 'quality')

###### FILTERING THE ROH DF
# which is probably a good idea, to filter the lower quality calls. 
filter_thresh <- 40
df <- df[df$quality >= filter_thresh, ]


# aggregate this to stats
df_stats_perchr <- data.frame(aggregate(df$length, by=list(df$sample, df$chr), FUN = sum),
                       roh_number = aggregate(df$length, by=list(df$sample, df$chr), FUN = length)$x)
df_stats_allchr <- data.frame(aggregate(df$length, by=list(df$sample), FUN = sum),
                       roh_number = aggregate(df$length, by=list(df$sample), FUN = length)$x)
colnames(df_stats_perchr) <- c('sample', 'chr', 'roh_total_length', 'roh_number')
colnames(df_stats_allchr) <- c('sample', 'roh_total_length', 'roh_number')

# TODO: use the basepairs of covered sequence instead of the total size
chr_sizes <- read.table(url('http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.chrom.sizes'))
colnames(chr_sizes) <- c('chr', 'size')
chr_sizes_vec <- setNames(chr_sizes$size, chr_sizes$chr)
total_chr_size <- sum(chr_sizes_vec[c('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22')])
df_stats_perchr$roh_total_length_norm <- 0
for (i in 1:nrow(df_stats_perchr)){
    df_stats_perchr[i, 'roh_total_length_norm'] <- df_stats_perchr[i, 'roh_total_length'] / chr_sizes_vec[paste0('chr', df_stats_perchr[i, 'chr'])]
}
df_stats_allchr$roh_total_length_norm <- df_stats_allchr$roh_total_length / total_chr_size
hist(df_stats_perchr$roh_total_length_norm)
hist(df_stats_allchr$roh_total_length_norm)
df_stats_perchr[order(df_stats_perchr$roh_total_length_norm, decreasing = T), ][1:20, ]

# annotate with some info from metadata
df_long <- df_stats_perchr
df_long <- df_long[df_long$sample != 'PD6567b',]
df_long$tumor_normal <- ega_metadata[df_long$sample, "tumor_normal"]
df_long$study <- ega_metadata[df_long$sample, "study"]
df_long$patient <- ega_metadata[df_long$sample, "patient"]
df_long$chr9p_loh <- ega_loh[df_long$patient, "chr9p_loh"]
df_long$any_loh <- ega_loh[df_long$patient, "any_loh"]
df_long$any_chr_event <- ega_loh[df_long$patient, "any_chr_event"]

# DISTRIBUTION OF CALL LENGTHS AND NUMBER
ggplot(df_long[df_long$chr==9,], aes(x=roh_number, y=roh_total_length, color=tumor_normal)) +
    geom_point() + 
    theme_bw() + 
    labs(title='ROH: EGA data, chr9, all samples')
ggplot(df_long[df_long$chr==9,], aes(x=roh_number, y=roh_total_length_norm, color=tumor_normal)) +
    geom_point() + 
    theme_bw() + 
    labs(title='ROH: EGA data, chr9, all samples')

p1 <- ggplot(df_long, aes(x=roh_number, y=roh_total_length, color=tumor_normal)) +
    geom_point(alpha=0.5) + 
    theme_bw() + 
    scale_y_log10() +
    facet_wrap(.~ chr) + 
    labs(title='ROH: EGA data, all samples')
p1
p2 <- ggplot(df_long, aes(x=roh_number, y=roh_total_length_norm)) +
    geom_point(alpha=0.25) + 
    theme_bw() + 
    facet_wrap(.~ chr) + 
    labs(title='ROH: EGA data, all samples')
p2

# BOXPLOTS COMPARED BETWEEN DISEASES: length
p1 <- ggplot(df_long, aes(x=study, y=roh_total_length, fill=study)) +
    geom_violin(trim=T)+
    geom_boxplot(width=0.1, fill="white")+
    theme_bw() + 
    scale_y_log10() +
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    labs(title='chr9 ROH, all samples', y="log10(roh_total_length)")
p2 <- ggplot(df_long[df_long$tumor_normal=='N' ,], aes(x=study, y=roh_total_length, fill=study)) +
    geom_violin(trim=T)+
    geom_boxplot(width=0.1, fill="white")+
    theme_bw() + 
    scale_y_log10() +
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    labs(title='chr9 ROH, germline samples', y="log10(roh_total_length)")
p3 <- ggplot(df_long[df_long$tumor_normal=='T',], aes(x=study, y=roh_total_length, fill=study)) +
    geom_violin(trim=T)+
    geom_boxplot(width=0.1, fill="white")+
    theme_bw() + 
    scale_y_log10() +
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    labs(title='chr9 ROH, tumor samples', y="log10(roh_total_length)")
plot_grid(p1,p2,p3, nrow=1)

# BOXPLOTS COMPARED BETWEEN DISEASES: number
p1 <- ggplot(df_long, aes(x=study, y=roh_number, fill=study)) +
    geom_violin(trim=T)+
    geom_boxplot(width=0.1, fill="white")+
    theme_bw() + 
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    labs(title='chr9 ROH, all samples')
p2 <- ggplot(df_long[df_long$tumor_normal=='N',], aes(x=study, y=roh_number, fill=study)) +
    geom_violin(trim=T)+
    geom_boxplot(width=0.1, fill="white")+
    theme_bw() + 
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    labs(title='chr9 ROH, germline samples')
p3 <- ggplot(df_long[df_long$tumor_normal=='T',], aes(x=study, y=roh_number, fill=study)) +
    geom_violin(trim=T)+
    geom_boxplot(width=0.1, fill="white")+
    theme_bw() + 
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV'))) + 
    labs(title='chr9 ROH, tumor samples')
plot_grid(p1,p2,p3, nrow=1)


######################################################################
# stratify by loh status
# all chromosomes
ggplot(df_long[!is.na(df_long$chr9p_loh) & 
               df_long$tumor_normal %in% c('N'),],
       aes(x=interaction(chr9p_loh, study), y=roh_total_length, fill=chr9p_loh)) +
    geom_boxplot() + 
    theme_bw() + 
    stat_compare_means(comparisons=list(c('FALSE.ET', 'TRUE.ET'),
                            c('FALSE.PMF', 'TRUE.PMF'),
                            c('FALSE.PV', 'TRUE.PV'))) + 
    facet_wrap(.~chr, ncol=6)+
    # ylim(c(0,2.75e8)) +
    labs(title='ROH: EGA data, per chromosome, germline samples')+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# just chr9
ggplot(df_long[!is.na(df_long$chr9p_loh) & 
               df_long$tumor_normal %in% c('N')& 
                   df_long$chr==9,],
       aes(x=interaction(chr9p_loh, study), y=roh_total_length, fill=chr9p_loh)) +
    geom_boxplot() + 
    theme_bw() + 
    stat_compare_means(comparisons=list(c('FALSE.ET', 'TRUE.ET'),
                            c('FALSE.PMF', 'TRUE.PMF'),
                            c('FALSE.PV', 'TRUE.PV'))) + 
    facet_wrap(.~chr, ncol=6)+
    # ylim(c(0,2.75e8)) +
    labs(title='ROH: EGA data, chr9, germline samples')+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# any chr event: No longer significant
ggplot(df_long[!is.na(df_long$any_chr_event) & 
               df_long$tumor_normal %in% c('N'),],
       aes(x=interaction(any_chr_event, study), y=roh_total_length, fill=any_chr_event)) +
    geom_boxplot() + 
    theme_bw() + 
    stat_compare_means(comparisons=list(c('FALSE.ET', 'TRUE.ET'),
                            c('FALSE.PMF', 'TRUE.PMF'),
                            c('FALSE.PV', 'TRUE.PV'))) + 
    facet_wrap(.~chr, ncol=6)+
    # ylim(c(0,2.75e8)) +
    labs(title='ROH: EGA data, per chromosome, germline samples') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


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
