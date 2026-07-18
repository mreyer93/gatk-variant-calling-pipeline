# Figure for the ET presentation about the changes in LOH over time
# based on data from the excel sheet of user

library(ggplot2)
library(cowplot)

loh_data <- read.table('~/local_data/LOH_reversal_figure/loh_reversal_data.tsv', sep='\t', quote='', header = T)
vaf_data <- read.table('~/local_data/LOH_reversal_figure/loh_reversal_VAF.tsv', sep='\t', quote='', header = T)

# change pre and post to be more specific
loh_data[loh_data$patient=='056-202' & loh_data$timepoint=='post', "timepoint"] <- 'Day 161'
loh_data[loh_data$patient=='030-201' & loh_data$timepoint=='post', "timepoint"] <- 'Day 171'
loh_data[loh_data$patient=='032-204' & loh_data$timepoint=='post', "timepoint"] <- 'Day 342'
loh_data[loh_data$patient=='003-202' & loh_data$timepoint=='post', "timepoint"] <- 'Day 166'
loh_data[loh_data$patient=='009-201' & loh_data$timepoint=='post', "timepoint"] <- 'Day 112'
loh_data[loh_data$timepoint=='pre', "timepoint"] <- 'Pre-treatment'

vaf_data[vaf_data$patient=='056-202' & vaf_data$timepoint=='post', "timepoint"] <- 'Day 161'
vaf_data[vaf_data$patient=='030-201' & vaf_data$timepoint=='post', "timepoint"] <- 'Day 171'
vaf_data[vaf_data$patient=='032-204' & vaf_data$timepoint=='post', "timepoint"] <- 'Day 342'
vaf_data[vaf_data$patient=='003-202' & vaf_data$timepoint=='post', "timepoint"] <- 'Day 166'
vaf_data[vaf_data$patient=='009-201' & vaf_data$timepoint=='post', "timepoint"] <- 'Day 112'
vaf_data[vaf_data$timepoint=='pre', "timepoint"] <- 'Pre-treatment'

loh_data$timepoint <- factor(loh_data$timepoint, levels=c('Pre-treatment', 'Day 161', 'Day 171', 'Day 112', 'Day 166', 'Day 342'))
vaf_data$timepoint <- factor(vaf_data$timepoint, levels=c('Pre-treatment', 'Day 161', 'Day 171', 'Day 112', 'Day 166', 'Day 342'))
loh_data$cell_type <- factor(loh_data$cell_type, levels=rev(c('wt', 'het', 'hom')))
pal <- c(wt='darkseagreen3', het='#E29585', hom='firebrick')

# change the names that appear in the facets
name_change <- c('003-202'='003-202\nJAK2V617F', 
                 '009-201'='009-201\nJAK2V617F',
                 '030-201'='030-201\nMPL_W515G',
                 '032-204'='032-204\nJAK2V617F',
                 '056-202'='056-202\nCALR_K385NCX')
loh_data$patient <- name_change[loh_data$patient]
vaf_data$patient <- name_change[vaf_data$patient]


p_celltype <- ggplot(loh_data, aes(x=timepoint, y=value/100, fill=cell_type)) + 
    geom_bar(position = 'stack', stat='identity') + 
    scale_fill_manual(values=pal, labels=c('Homozygous mutant','Heterozygous mutant' ,'Wild-type')) +
    facet_wrap(.~patient, nrow=1, scales='free_x') + 
    labs(x='Timepoint', y='Percentage') +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    theme_bw() + 
    guides(fill=guide_legend(title='Granulocyte\ngenotype', reverse = T))

p_vaf <- ggplot(vaf_data, aes(x=timepoint, y=vaf/100)) + 
    geom_bar(stat='identity', fill='grey20') +
    facet_wrap(.~patient, nrow=1, scales='free_x') + 
    labs(title='Granulocyte genotypes: pre- and post-treatment', y='Mutant Allele Frequency') +
    theme_bw() + 
    scale_x_discrete(drop=T) +
    scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank()) 


p3 <- plot_grid(
    plot_grid(
        p_vaf
        , NULL
        , p_celltype + theme(legend.position = "none")
        , ncol = 1
        , align = "hv",
        rel_heights = c(1,-0.15,1))
    , plot_grid(
        ggplot() + theme_minimal(),
        get_legend(p_celltype)
        , ncol =1)
    , rel_widths = c(6,1.2)
)
p3

# save a combined plot
pdf('~/local_data/LOH_reversal_figure/LOH_reversal_figure_FINAL.pdf', height=6, width = 11)
print(p3)
dev.off()

png('~/local_data/LOH_reversal_figure/LOH_reversal_figure_FINAL.png', height=6, width = 11, res=300, units="in")
print(p3)
dev.off()


