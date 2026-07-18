library(ggplot2)
library(ggpubr)
library(cowplot)

# ega metadata joined with the ROH results
df <- read.table('~/local_data/ega_metadata_join_roh.tsv', sep='\t', quote='', header=T, comment.char = '')

p1 <- ggplot(df, aes(x=roh_number, y=roh_total_length, color=tumor_normal)) + 
    geom_point() + 
    theme_bw() + 
    labs(title='EGA chr9 ROH stats')
p2 <- ggplot(df, aes(x=roh_number, y=roh_total_length, color=study)) + 
    geom_point() + 
    theme_bw() +
    labs(title='')
plot_grid(p1, p2, ncol = 2)

p1 <- ggplot(df, aes(x=tumor_normal, y=roh_total_length, fill=tumor_normal)) + 
    geom_boxplot() + 
    theme_bw() +
    stat_compare_means() + 
    labs(title='EGA chr9 ROH stats')
p2 <- ggplot(df, aes(x=tumor_normal, y=roh_number, fill=tumor_normal)) + 
    geom_boxplot() + 
    theme_bw() +
    stat_compare_means() +
    labs(title='')
plot_grid(p1,p2)


p1 <- ggplot(df, aes(x=study, y=roh_total_length, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='EGA chr9 ROH stats') +
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
p2 <- ggplot(df, aes(x=study, y=roh_number, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='')+
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
plot_grid(p1,p2)

p1 <- ggplot(df[df$tumor_normal=='N',], aes(x=study, y=roh_total_length, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='EGA chr9 ROH stats: NORMAL SAMPLES') +
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
p2 <- ggplot(df[df$tumor_normal=='N',], aes(x=study, y=roh_number, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='')+
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
p3 <- ggplot(df[df$tumor_normal=='T',], aes(x=study, y=roh_total_length, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='EGA chr9 ROH stats: TUMOR SAMPLES') +
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
p4 <- ggplot(df[df$tumor_normal=='T',], aes(x=study, y=roh_number, fill=study)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='')+
    stat_compare_means(comparisons = list(c('ET', 'PMF'), 
                                          c('PMF', 'PV'),
                                          c('ET', 'PV')))
plot_grid(p1,p2,p3,p4, ncol=2)

summary(lm(roh_number~tumor_normal, data=df))
summary(lm(roh_total_length~tumor_normal, data=df))
summary(lm(roh_number~study, data=df))
summary(lm(roh_total_length~study, data=df))

# what if we join with the LOH metadata from the tumor samples
tumor_metadata <- read.table('~/local_data/ega_loh_metadata.tsv', sep='\t', quote='' ,header = T)
# annotate other df with this
library(dplyr)
df_join <- left_join(df, tumor_metadata, by='patient')
# this will only be of the paired samples, the other data will be NA


p1 <- ggplot(df_join[!is.na(df_join$chr9p_loh) & 
                         df_join$tumor_normal %in% c('N'),],
             aes(x=interaction(any_chr_event, study), y=roh_total_length, fill=any_chr_event)) + 
    geom_boxplot() + 
    theme_bw() +
    labs(title='EGA chr9 ROH: normal samples, annotated with tumor LOH status') +
    stat_compare_means(comparisons = list(c('FALSE.ET', 'TRUE.ET'),
                                          c('FALSE.PMF', 'TRUE.PMF'),
                                          c('FALSE.PV', 'TRUE.PV')))
p1


# what about the "p6" data? 
p6 <- read.table('~/local_data/EGA_haplotypecaller_p6_calculated_integers.tsv', sep='\t', quote='', header=T, check.names = F)
rownames(df_join) <- df_join$sample
df_join$p6_num_hom <- colSums(p6[,df_join$sample] ==2 )
df_join$p6_num_hom_5 <- df_join$p6_num_hom==5

p1 <- ggplot(df_join, aes(x=roh_number, y=roh_total_length, color=p6_num_hom_5)) + 
    geom_point() + 
    theme_bw() + 
    facet_wrap(.~tumor_normal)+ 
    labs(title='EGA chr9 ROH stats')
p1


df_join[, c('patient', 'tumor_normal', "p6_num_hom")]
table(df_join[, c('tumor_normal', "p6_num_hom")])
table(df_join[, c('study', 'tumor_normal', "p6_num_hom")])
table(df_join[, c('study', 'tumor_normal' )])


p1 <- ggplot(df_join, aes(x=p6_num_hom_5, y=roh_total_length, fill=p6_num_hom_5)) + 
    geom_boxplot() + 
    theme_bw() +
    stat_compare_means() + 
    facet_wrap(~study) +
    labs(title='EGA chr9 ROH stats')
p1

p2 <- ggplot(df_join, aes(x=tumor_normal, y=roh_number, fill=tumor_normal)) + 
    geom_boxplot() + 
    theme_bw() +
    stat_compare_means() +
    labs(title='')
plot_grid(p1,p2)