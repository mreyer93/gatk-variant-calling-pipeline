# investigating ROH in EGA data using bcftools ROH 
# with the RG data instead of the ST tags
# should be better? 
library(ggplot2)
library(ggExtra)
library(ggpubr)
library(cowplot)
# get the EGA sample annotations from the metadata file
ega_metadata <- read.table('~/local_data/EGA_data/ega_metadata_join_roh.tsv', sep='\t', header = T, comment.char = '', quote='')
rownames(ega_metadata) <- ega_metadata$sample
ega_loh <- read.table('~/local_data/EGA_data/ega_loh_metadata.tsv', sep='\t', header = T, comment.char = '', quote='')
rownames(ega_loh) <- ega_loh$patient
ega_loh$chr9_loh <- ega_loh$chr9p_loh
ega_loh["PDNEW16","chr9_loh"] <- TRUE

sample_list <-  read.table('/path/to/data/p6/EGA_haplotypecaller/sample_list.txt')[,1]
sample_list <- sample_list[sample_list != 'PDNEW36a']
sample_list <- sample_list[sample_list != 'PD6567b']
basedir <- '/path/to/data/p6/EGA_haplotypecaller/'

roh_stats_df_ST <- read.table('/path/to/data/p6/EGA_haplotypecaller/roh_stats_df/chr9_roh_stats.txt', sep='\t', quote='', header=T, row.names = 1)
colnames(roh_stats_df_ST)[1] <- 'sample'
head(roh_stats_df_ST[order(roh_stats_df_ST$roh_total_length, decreasing = T), ])

outdir <- '/path/to/data/p6/EGA_haplotypecaller/roh_stats_df_RG/'
# for (chr in 1:22){
for (chr in 9){
    df_chr <- data.frame()
    print(chr)    
    for (sample in sample_list){
        print(sample)
        gzf <- gzfile(file.path(basedir, paste0("chr", chr), paste0(sample, "_roh.txt.gz")))
        lines <- readLines(gzf)
        lines_ind <- sapply(lines, function(x) substr(x,1,2)=='RG')
        lines_subset <- lines[lines_ind]
        df <- read.table(sep='\t', header=F, text=lines_subset)
        colnames(df) <- c('class', 'sample', 'chr', 'start', 'end', 'length', 'num_markers', 'quality')
        df_chr <- rbind(df_chr, df)
    }
    write.table(df_chr, file.path(outdir, paste0('chr', chr, '_roh_stats.txt')), sep='\t', quote=F, row.names=F, col.names=T)

}
df_chr$patient <- ega_metadata[df_chr$sample, "patient"]
df_chr$tumor_normal <- ega_metadata[df_chr$sample, "tumor_normal"]
df_chr$chr9p_loh <- ega_loh[df_chr$patient, "chr9p_loh"]
df_chr$chr9_loh <- ega_loh[df_chr$patient, "chr9_loh"]


p <- ggplot(df_chr, aes(x=quality, y=length)) +
    geom_point(alpha=0.25) + 
    scale_y_log10()+
    labs(title='EGA chr9, all ROH events')
ggMarginal(p, type='histogram')

p <- ggplot(df_chr, aes(x=quality, y=length)) +
    geom_point(alpha=0.25) + 
    labs(title='EGA chr9, all ROH events')
ggMarginal(p, type='histogram')


# filter on some intelligent and somewhat random number, like 40
df_chr_filt <- df_chr[df_chr$quality >=40, ]

View(df_chr_filt[order(df_chr_filt$length, decreasing = T), ])
p <- ggplot(df_chr_filt, aes(x=quality, y=length, color=chr9_loh, shape=tumor_normal)) +
    geom_point() + 
    scale_y_log10() + 
    theme_bw() 
ggMarginal(p, type='density')

ggplot(df_chr_filt, aes(x=start, y=length, color=chr9_loh, shape=tumor_normal)) +
    geom_point(size=2.5) + 
    theme_bw()


# aggregate by length
sample_sum <- aggregate(df_chr_filt$length, by=list(df_chr_filt$sample), FUN=sum)
sample_num <- aggregate(df_chr_filt$length, by=list(df_chr_filt$sample), FUN=length)
agg_df <- data.frame(sample=sample_sum$Group.1, roh_total_length=sample_sum$x, roh_number=sample_num$x)
agg_df$patient <- ega_metadata[agg_df$sample, "patient"]
agg_df$tumor_normal <- ega_metadata[agg_df$sample, "tumor_normal"]
agg_df$study <- ega_metadata[agg_df$sample, "study"]
agg_df$chr9p_loh <- ega_loh[agg_df$patient, "chr9p_loh"]
agg_df$any_loh <- ega_loh[agg_df$patient, "any_loh"]
agg_df$any_chr_event <- ega_loh[agg_df$patient, "any_chr_event"]
rownames(agg_df) <- agg_df$sample

ggplot(agg_df, aes(x=roh_number, y=roh_total_length, color=tumor_normal)) + 
    geom_point() + 
    theme_bw() + 
    labs(title='EGA data filtered ROH, chr9')

p1 <- ggplot(agg_df, aes(x=tumor_normal, y=roh_total_length, fill=tumor_normal)) + 
    geom_boxplot() + 
    theme_bw() +
    stat_compare_means() + 
    labs(title='EGA chr9 ROH stats')
p2 <- ggplot(agg_df, aes(x=tumor_normal, y=roh_number, fill=tumor_normal)) + 
    geom_boxplot() + 
    theme_bw() +
    stat_compare_means() +
    labs(title='')
plot_grid(p1,p2)

# difference between study
p1 <- ggplot(agg_df, aes(x=study, y=roh_total_length, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='EGA chr9 ROH stats') +
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
p2 <- ggplot(agg_df, aes(x=study, y=roh_number, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='')+
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
plot_grid(p1,p2)


# normal samples only
p1 <- ggplot(agg_df[agg_df$tumor_normal=='N', ], aes(x=study, y=roh_total_length, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='EGA chr9 ROH stats') +
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
p2 <- ggplot(agg_df[agg_df$tumor_normal=='N', ], aes(x=study, y=roh_number, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='')+
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
plot_grid(p1,p2)


# association with LOH
p1 <- ggplot(agg_df[!is.na(agg_df$chr9p_loh) & 
                        agg_df$tumor_normal %in% c('N'),],
             aes(x=interaction(any_chr_event, study), y=roh_total_length, fill=any_chr_event)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='EGA chr9 ROH: normal samples, annotated with tumor LOH status') +
    stat_compare_means(comparisons = list(c('FALSE.ET', 'TRUE.ET'),
                                          c('FALSE.PMF', 'TRUE.PMF'),
                                          c('FALSE.PV', 'TRUE.PV')))
p1

# what about this other method? 
adf <- data.frame(sample=sample_list, audacity_sum = 0)
rownames(adf) <- adf$sample
for (sample in sample_list){
    print(sample)
    sample_f <- paste0('~/data/p6/EGA_haplotypecaller/audacity/chr9/', sample, '/', sample, '_DIDOH3M2Regions.txt')
    sample_df <- read.table(sample_f, sep='\t', quote='', header=F)
    adf[sample, "audacity_sum"] <- sum(sample_df$V4)    
}

hist(adf$audacity_sum)

adf$bcftools_sum <- agg_df[adf$sample, "roh_total_length"]
ggplot(adf, aes(x=bcftools_sum, y=audacity_sum)) + 
    geom_point()
summary(lm(audacity_sum~bcftools_sum, data=adf))
