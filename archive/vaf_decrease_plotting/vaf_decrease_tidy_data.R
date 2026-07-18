# read in the table from the VAF data crunching
library(ggplot2)
library(reshape2)
library(cowplot)
library(ggrepel)
library(dplyr)

vaf <- read.table('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/tidy_vaf_data.tsv',
                  sep='\t', quote='', header=T)
loh <- read.table('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/loh_events_patient.tsv',
                  sep='\t', quote='', header=T, comment.char = '')
loh.jak2.pts <- loh[loh$X9p_loh,"patient"]

# some data cleaning
# remove patients without DX sample
pt.dx <- unique(vaf[vaf$timepoint=='Dx', "patient"])
vaf <- vaf[vaf$patient %in% pt.dx, ]

# with at least one time course sample
vaf$patient_timepoint <- paste0(vaf$patient, '_', vaf$timepoint)
vaf.sub <- vaf[!duplicated(vaf$patient_timepoint), ]
keep.patients <- names(which(table(vaf.sub$patient) >1))
vaf <- vaf[vaf$patient %in% keep.patients,]
rownames(vaf)<- 1:nrow(vaf)
# only keep tumor
vaf <- vaf[vaf$sample_type=="T",]
# remove asxl1 642 variants because they're likely false positives
vaf <- vaf[vaf$mut!='ASXL1_-642X', ]
# remove duplicate MPL variants
vaf <- vaf[vaf$mut!='MPL_W515R', ]

# aggregate to mean value per tissue
vaf$patient_timepoint_mutation <- paste(vaf$patient, vaf$timepoint, vaf$mut, sep='_')
agg.vaf <- aggregate(vaf$VAF, list(vaf$patient_timepoint_mutation), mean)
vaf.agg <- vaf[!duplicated(vaf$patient_timepoint_mutation), ]
rownames(agg.vaf) <- agg.vaf$Group.1
vaf.agg$VAF <- agg.vaf[vaf.agg$patient_timepoint_mutation,"x"]
vaf.agg <- vaf.agg[, colnames(vaf.agg) != 'tissue']
vaf <- vaf.agg
vaf$patient_mutation <- paste(vaf$patient, vaf$mut, sep='_')

# done with CLEANING
sort(table(vaf$mut), decreasing = T)

plot.mut <- 'JAK2_V617F'
plot.mut <- 'CALR_52b_del'
plot.mut <- 'ASXL1_R693*'
plot.mut <- 'MPL_W515*'
plot.mut <- 'MPL_W515R'
ggplot(vaf[vaf$mut==plot.mut & vaf$sample_type=="T",], aes(x=days, y=VAF, group=patient)) +
    geom_line() +
    labs(title=plot.mut, x='Days of treatment', y='Variant allele frequency') +
    ylim(c(0,100)) +
    theme_bw()



# do relative calculations
vaf.initial <- rep(0, nrow(vaf))
for (i in 1:nrow(vaf)){
    print(i)
    this.row <- vaf[i, ]
    if(this.row$days==0){
        vaf.initial[i] <- 0
    } else {
        match.row <- vaf[vaf$patient==this.row$patient &
                         # vaf$tissue==this.row$tissue &
                         vaf$mut==this.row$mut &
                         vaf$days==0,]
        if(nrow(match.row) ==0){
            print(match.row)
            vaf.initial[i] <- NA
        } else {
            vaf.initial[i] <- match.row$VAF
        }
    }
}

# calculate relative changes of these
    vaf$vaf.initial <- vaf.initial
vaf$pct.change <- (vaf$VAF - vaf$vaf.initial) / vaf$vaf.initial * 100
vaf$pct.change[vaf$vaf.initial==0] <- 0
vaf$fold.change <- vaf$VAF / vaf$vaf.initial
vaf$fold.change[vaf$vaf.initial==0] <- 1
vaf$absolute.change <- vaf$VAF - vaf$vaf.initial
vaf$absolute.change[vaf$vaf.initial==0] <- 0
vaf <- vaf[!is.na(vaf$pct.change),]

# metric of increase/decrease
major.decrease <- -25
major.increase <- 25
stable <- 5
# now doing just 10% thresholds
increase <- 10
decrease <- -10
# value for last sample from a patient
vaf.final <- vaf[order(vaf$days, decreasing = T), ]
vaf.final <- vaf.final[!duplicated(vaf.final$patient_mutation),]
vaf.final$classification <- 'minor decrease'
vaf.final[vaf.final$pct.change < major.decrease, "classification"] <- 'major decrease'
vaf.final[vaf.final$pct.change > stable & vaf.final$pct.change <= major.increase, "classification"] <- 'minor increase'
vaf.final[vaf.final$pct.change < stable & vaf.final$pct.change >-stable, "classification"] <- 'stable'
vaf.final[vaf.final$pct.change < major.increase & vaf.final$pct.change >= stable, "classification"] <- 'minor increase'
vaf.final[vaf.final$pct.change > major.increase, "classification"] <- 'major increase'

vaf.final$classification.simple <- 'Stable'
vaf.final[vaf.final$pct.change <= decrease, "classification.simple"] <- 'Decrease'
vaf.final[vaf.final$pct.change >= increase, "classification.simple"] <- 'Increase'


patient.mut.classification <- vaf.final$classification
names(patient.mut.classification) <- vaf.final$patient_mutation
vaf$classification <- patient.mut.classification[vaf$patient_mutation]

patient.mut.classification.simple <- vaf.final$classification.simple
names(patient.mut.classification.simple) <- vaf.final$patient_mutation
vaf$classification.simple <- patient.mut.classification.simple[vaf$patient_mutation]

vaf$label_final <- ''
vaf[vaf$patient_timepoint_mutation %in% vaf.final$patient_timepoint_mutation, 'label_final'] <- vaf[vaf$patient_timepoint_mutation %in% vaf.final$patient_timepoint_mutation, "patient"]

pal = c( 'major decrease'='steelblue4', 'minor decrease'='lightskyblue3', stable='grey50', 'minor increase' = 'lightcoral', 'major increase' = 'firebrick')
pal.simple = c( 'Decrease'='steelblue4', 'Stable'='grey50', 'Increase' = 'firebrick')

# do this for all mutations and export to PDF
pdf('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/CTP-102_vaf_over_time_mutation_relative.pdf', height=5, width=7)
plot.muts <- sort(table(vaf.final$mut), decreasing = T)
for (plot.mut in names(plot.muts)[plot.muts>1]){
plotme <- vaf[vaf$mut==plot.mut & vaf$sample_type=="T",]
plotme[plotme$days >300, "days"] <- 300
plotme[plotme$pct.change >100, "pct.change"] <- 100
# if(plot.mut =='JAK2_V617F'){
#     plotme$GC <- plotme$patient %in% loh.jak2.pts
#     p1 <- ggplot(plotme, aes(x=days, y=pct.change, col=classification, group=patient, label=label_final, lty=GC)) +
#         geom_line(size=0.9) +
#         scale_color_manual(values=pal) +
#         geom_hline(yintercept = 0, lty=2, col='grey80') +
#         labs(title=paste('Precent change: ', plot.mut), x='Days of treatment', y='Variant allele frequency') +
#         # geom_text_repel(box.padding = 0.5, max.overlaps = 10) +
#         theme_bw()
# }
p1 <- ggplot(plotme, aes(x=days, y=pct.change/100, col=classification, group=patient, label=label_final)) +
    geom_line(size=0.9) +
    scale_color_manual(values=pal) +
    geom_hline(yintercept = 0, lty=2, col='grey80') +
    labs(title=paste('Variant allele frequency (VAF):', plot.mut), subtitle='Mean value per timepoint, relative to enrollment',
         x='Days of treatment', y='VAF % change relative to enrollment', color='Classification') +
    scale_y_continuous(labels=scales::percent)  +
    # geom_text_repel(box.padding = 0.5, max.overlaps = 10) +
    theme_bw()
p1.simple <- ggplot(plotme, aes(x=days, y=pct.change/100, col=classification.simple, group=patient, label=label_final)) +
    geom_line(size=0.9) +
    scale_color_manual(values=pal.simple) +
    geom_hline(yintercept = 0, lty=2, col='grey80') +
    labs(title=paste('Variant allele frequency (VAF):', plot.mut), subtitle='Mean value per timepoint, relative to enrollment',
         x='Days of treatment', y='VAF % change relative to enrollment', color='Classification') +
    scale_y_continuous(labels=scales::percent)  +
    # geom_text_repel(box.padding = 0.5, max.overlaps = 10) +
    theme_bw()
print(p1)
print(p1.simple)
}
dev.off()

# absolute percentages
pdf('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/CTP-102_vaf_over_time_mutation_absolute.pdf', height=5, width=7)
plot.muts <- sort(table(vaf.final$mut), decreasing = T)
for (plot.mut in names(plot.muts)[plot.muts>1]){
plotme <- vaf[vaf$mut==plot.mut & vaf$sample_type=="T",]
plotme[plotme$days >300, "days"] <- 300
plotme[plotme$pct.change >100, "pct.change"] <- 100

p1 <- ggplot(plotme, aes(x=days, y=VAF, group=patient, color=classification, label=label_final)) +
    geom_line(size=0.9) +
    scale_color_manual(values=pal) +
    labs(title=plot.mut, x='Days of treatment', y='Variant allele frequency', color='Classification') +
    ylim(c(0,100)) +
    theme_bw()
p1.simple <- ggplot(plotme, aes(x=days, y=VAF, group=patient, color=classification.simple, label=label_final)) +
    geom_line(size=0.9) +
    scale_color_manual(values=pal.simple) +
    labs(title=plot.mut, x='Days of treatment', y='Variant allele frequency', color='Classification') +
    ylim(c(0,100)) +
    theme_bw()
print(p1)
print(p1.simple)
}
dev.off()

# other stats you can plot
ggplot(plotme, aes(x=days, y=log2(fold.change), group=patient)) +
    geom_line() +
    labs(title=paste('log2 fold change: ', plot.mut), x='Days of treatment', y='Variant allele frequency') +
    theme_bw()

ggplot(plotme, aes(x=days, y=absolute.change, group=patient)) +
    geom_line() +
    labs(title=paste('percentage point change ', plot.mut), x='Days of treatment', y='Variant allele frequency') +
    theme_bw()
ggplot(plotme, aes(x=days, y=VAF, group=patient, color=classification, label=label_final, lty=GC)) +
    geom_line() +
    scale_color_manual(values=pal) +
    labs(title=plot.mut, x='Days of treatment', y='Variant allele frequency') +
    ylim(c(0,100)) +
    theme_bw()

### patient centric view with multiple mutations tracked over time?
# do this for all patients
pdf('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/CTP-102_vaf_over_time_patient.pdf')
for (plot.patient in unique(vaf$patient)){
    print(plot.patient)
    plotme <- vaf[vaf$patient == plot.patient & vaf$sample_type=="T",]
    plotme[plotme$days >300, "days"] <- 300
    plotme[plotme$pct.change >100, "pct.change"] <- 100
    # absolute VAF
    p1 <- ggplot(plotme, aes(x=days, y=VAF, color=classification.simple, group=patient)) +
        geom_line(size=0.9) +
        scale_color_manual(values=pal.simple) +
        facet_wrap(.~mut) +
        labs(title=paste('Patient', plot.patient), x='Days of treatment', y='Variant allele frequency') +
        theme_bw() +
        ylim(c(0,100))

    # relative VAF
    p2 <- ggplot(plotme, aes(x=days, y=pct.change/100, color=classification.simple, group=patient)) +
        geom_line(size=0.9) +
        scale_color_manual(values=pal.simple) +
        facet_wrap(.~mut) +
        labs(
             x='Days of treatment',
             y='VAF % change from baseline') +
        scale_y_continuous(labels=scales::percent)  +
        theme_bw()

    print(plot_grid(p1, p2, nrow=2))
}
dev.off()

#  get stats behind each of these mutations, number of patients in each catergory, etc.
# total patients with time course tumor samples
pdf('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/CTP-102_patient_classification_fraction_bar.pdf', height=3, width=8)
pielist <- list()
for (plot.mut in names(plot.muts)[plot.muts>1]){
    jak2.patients <- vaf[!duplicated(vaf$patient_mutation) & vaf$mut==plot.mut, c('patient', 'classification'), ]
    jak2.patients <- jak2.patients[order(jak2.patients$classification), ]
    jak2.table <- as.data.frame(sort(table(vaf[!duplicated(vaf$patient_mutation) & vaf$mut==plot.mut, 'classification' ]), decreasing = T))
    jak2.table$mut <- plot.mut
    colnames(jak2.table)[1] <- 'classification'
    rownames(jak2.table) <- jak2.table$classification
    jak2.table <- jak2.table[c('major decrease', 'minor decrease', 'stable', 'minor increase', 'major increase'),]
    jak2.table$classification <- factor(jak2.table$classification, levels=rev(c('major decrease', 'minor decrease', 'stable', 'minor increase', 'major increase')))
    jak2.table$pos <- cumsum(jak2.table$Freq) - (jak2.table$Freq/2)

    jak2.table.simple <- as.data.frame(sort(table(vaf[!duplicated(vaf$patient_mutation) & vaf$mut==plot.mut, 'classification.simple' ]), decreasing = T))
    jak2.table.simple$mut <- plot.mut
    colnames(jak2.table.simple)[1] <- 'classification.simple'
    rownames(jak2.table.simple) <- jak2.table.simple$classification.simple
    jak2.table.simple <- jak2.table.simple[c('Decrease', 'Stable', 'Increase'),]
    jak2.table.simple$classification.simple <- factor(jak2.table.simple$classification.simple, levels=rev(c('Decrease', 'Stable', 'Increase')))
    jak2.table.simple$pos <- cumsum(jak2.table.simple$Freq) - (jak2.table.simple$Freq/2)

    p <- ggplot(jak2.table, aes(x=mut, y=Freq, fill=classification)) +
        geom_bar(stat='identity', position='stack') +
        scale_fill_manual(values=pal) +
        theme_bw() +
        coord_flip() +
        labs(title=paste(plot.mut, "patient classification"), x='', y='', fill='Classification') +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        geom_text(aes(y=pos, label=Freq), nudge_x = 0.52, size=5)
    p.simple <- ggplot(jak2.table.simple, aes(x=mut, y=Freq, fill=classification.simple)) +
        geom_bar(stat='identity', position='stack') +
        scale_fill_manual(values=pal.simple) +
        theme_bw() +
        coord_flip() +
        labs(title=paste(plot.mut, "patient classification"), x='', y='', fill='Classification') +
        theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank()) +
        geom_text(aes(y=pos, label=Freq), nudge_x = 0.52, size=5)
    p.pie <- p + coord_polar('y', start=0) + theme_void()
    p.simple.pie <- p.simple + coord_polar('y', start=0) + theme_void()
    print(p)
    print(p.simple)
    pielist <- append(pielist, list(p.pie))
    pielist <- append(pielist, list(p.simple.pie))
}
dev.off()

pdf('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/CTP-102_patient_classification_fraction_pie.pdf', height=4, width=6)
for (i in 1:length(pielist)){
    print(pielist[[i]])
}
dev.off()

# looking for correlations with clinical data
# read in the sheet which was exported from excel
clin <- read.table('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/ctp_102_patient_data_export.tsv', sep='\t', quote='', header=T)

# for each patient, calculate some stats from inital to final
# symptom score
# spleen volume
# spleen size
# fatigue
# parsed all this data in a spreadsheet and now have the pct changes
clin.change <- read.table('~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/CTP-102_clinical_data_pct_change.tsv', sep='\t', quote='', header=T)

# can I do a straight regression?
plot.mut <- 'JAK2_V617F'
plot.mut <- 'CALR_52b_del'
vaf.final.join <- left_join(vaf.final, clin.change, by = 'patient')
write.table(vaf.final.join, '~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/CTP-102_vaf_final_join_clinical.tsv', sep='\t', quote=F, col.names=T, row.names=F)

vj2 <- vaf.final[vaf.final$mut == plot.mut, ]
vj2.clin <- left_join(vj2, clin.change, by = 'patient')
# what if I exclude that outlier 200%
vj2.clin <- vj2.clin[vj2.clin$pct.change <150, ]
p1 <- ggplot(vj2.clin, aes(x=pct.change, y=spleen_vol_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p2 <- ggplot(vj2.clin, aes(x=vaf.initial, y=spleen_vol_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p3 <- ggplot(vj2.clin, aes(x=absolute.change, y=spleen_vol_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p4 <- ggplot(vj2.clin, aes(x=pct.change, y=spleen_len_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p5 <- ggplot(vj2.clin, aes(x=vaf.initial, y=spleen_len_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p6 <- ggplot(vj2.clin, aes(x=absolute.change, y=spleen_len_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p7 <- ggplot(vj2.clin, aes(x=pct.change, y=mpn10_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p8 <- ggplot(vj2.clin, aes(x=vaf.initial, y=mpn10_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p9 <- ggplot(vj2.clin, aes(x=absolute.change, y=mpn10_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p10 <- ggplot(vj2.clin, aes(x=pct.change, y=fatigue_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p11 <- ggplot(vj2.clin, aes(x=vaf.initial, y=fatigue_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p12 <- ggplot(vj2.clin, aes(x=absolute.change, y=fatigue_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()

p13 <- ggplot(vj2.clin, aes(x=pct.change, y=WBC_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p14 <- ggplot(vj2.clin, aes(x=vaf.initial, y=WBC_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p15 <- ggplot(vj2.clin, aes(x=absolute.change, y=WBC_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()

p16 <- ggplot(vj2.clin, aes(x=pct.change, y=neutro_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p17 <- ggplot(vj2.clin, aes(x=vaf.initial, y=neutro_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()
p18 <- ggplot(vj2.clin, aes(x=absolute.change, y=neutro_change)) + geom_point() + geom_smooth(method='lm', se=T) + theme_bw()

plotlist <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18)
title <- ggdraw() + draw_label(paste('Clinical correlations: ', plot.mut), fontface='bold')

plot_grid(title, plot_grid(plotlist = plotlist, nrow=6), ncol=1, rel_heights=c(0.06, 1)) # rel_heights values control title margins


summary(lm(WBC_change ~ pct.change, data=vj2.clin))
summary(lm(WBC_change ~ vaf.initial, data=vj2.clin))
summary(lm(WBC_change ~ absolute.change, data=vj2.clin))
summary(lm(WBC_change ~ classification, data=vj2.clin))
summary(lm(WBC_change ~ classification.simple, data=vj2.clin))

summary(lm(neutro_change ~ pct.change, data=vj2.clin))
summary(lm(neutro_change ~ vaf.initial, data=vj2.clin))
summary(lm(neutro_change ~ absolute.change, data=vj2.clin))
summary(lm(neutro_change ~ classification, data=vj2.clin))
summary(lm(neutro_change ~ classification.simple, data=vj2.clin))

summary(lm(spleen_vol_change ~ pct.change, data=vj2.clin))
summary(lm(spleen_vol_change ~ vaf.initial, data=vj2.clin))
summary(lm(spleen_vol_change ~ absolute.change, data=vj2.clin))
summary(lm(spleen_vol_change ~ classification, data=vj2.clin))
summary(lm(spleen_vol_change ~ classification.simple, data=vj2.clin))

summary(lm(spleen_len_change ~ pct.change, data=vj2.clin))
summary(lm(spleen_len_change ~ vaf.initial, data=vj2.clin))
summary(lm(spleen_len_change ~ absolute.change, data=vj2.clin))
summary(lm(spleen_len_change ~ classification, data=vj2.clin))
summary(lm(spleen_len_change ~ classification.simple, data=vj2.clin))

summary(lm(mpn10_change ~ pct.change, data=vj2.clin))
summary(lm(mpn10_change ~ vaf.initial, data=vj2.clin))
summary(lm(mpn10_change ~ absolute.change, data=vj2.clin))
summary(lm(mpn10_change ~ classification, data=vj2.clin))
summary(lm(mpn10_change ~ classification.simple, data=vj2.clin))

summary(lm(fatigue_change ~ pct.change, data=vj2.clin))
summary(lm(fatigue_change ~ vaf.initial, data=vj2.clin))
summary(lm(fatigue_change ~ absolute.change, data=vj2.clin))
summary(lm(fatigue_change ~ classification, data=vj2.clin))
summary(lm(fatigue_change ~ classification.simple, data=vj2.clin))

classification.map = c("stable" =0, "minor decrease" = -1,  "major decrease" = -2, "major increase" =2, "minor increase" =1)
cor(vj2.clin$spleen_vol_change, vj2.clin$pct.change, method='spear', use='complete.obs')
cor(vj2.clin$spleen_vol_change, vj2.clin$vaf.initial, method='spear', use='complete.obs')
cor(vj2.clin$spleen_vol_change, vj2.clin$absolute.change, method='spear', use='complete.obs')
cor(vj2.clin$spleen_vol_change, classification.map[vj2.clin$classification], method='spear', use='complete.obs')
cor(vj2.clin$spleen_len_change, vj2.clin$pct.change, method='spear', use='complete.obs')
cor(vj2.clin$spleen_len_change, vj2.clin$vaf.initial, method='spear', use='complete.obs')
cor(vj2.clin$spleen_len_change, vj2.clin$absolute.change, method='spear', use='complete.obs')
cor(vj2.clin$spleen_len_change, classification.map[vj2.clin$classification], method='spear', use='complete.obs')
cor(vj2.clin$mpn10_change, vj2.clin$pct.change, method='spear', use='complete.obs')
cor(vj2.clin$mpn10_change, vj2.clin$vaf.initial, method='spear', use='complete.obs')
cor(vj2.clin$mpn10_change, vj2.clin$absolute.change, method='spear', use='complete.obs')
cor(vj2.clin$mpn10_change, classification.map[vj2.clin$classification], method='spear', use='complete.obs')
cor(vj2.clin$fatigue_change, vj2.clin$pct.change, method='spear', use='complete.obs')
cor(vj2.clin$fatigue_change, vj2.clin$vaf.initial, method='spear', use='complete.obs')
cor(vj2.clin$fatigue_change, vj2.clin$absolute.change, method='spear', use='complete.obs')
cor(vj2.clin$fatigue_change, classification.map[vj2.clin$classification], method='spear', use='complete.obs')

# write out date from final timepoint for all mutations
vaf.export <- vaf.final[vaf.final$mut %in% names(plot.muts)[plot.muts>1], ]
vaf.export <- vaf.export[order(vaf.export$pct.change),]
write.table(vaf.export, '~/pcloud_sync/project/vaf_decrease_plotting/second_maf_2021-11-06/vaf_final.tsv', sep = '\t', quote=F, col.names = T, row.names = F)

# Todo
# repeat figures with 10% thresholds
# make bar/pie charts with current and 10% threshold
# send vaf final with clinical data smashed in.
