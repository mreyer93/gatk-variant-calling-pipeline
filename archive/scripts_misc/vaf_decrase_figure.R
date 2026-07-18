# plotting VAF decrease over time
# in a single patient vs multiple patients
library(ggplot2)
library(reshape2)

# JAK2 V617F ONLY
vaf.jak2 <- read.table("~/pcloud_sync/project/vaf_decrease_plotting/jak2v617f_data.tsv", sep='\t', quote='', header=T)
# all variants in long form
vaf.all <- read.table("~/pcloud_sync/project/vaf_decrease_plotting/all_variant_data.tsv", sep='\t', quote='', header=T)

################# SOME QC THINGS #################
# limit to patients with multiple timepoints
vaf.all <- vaf.all[vaf.all$fo_tumor >0, ]
# remove asxl1 642 variants because they're likely false positives
vaf.all <- vaf.all[vaf.all$variant!='ASXL1_-642X', ]

# how many patients per variant?
vaf.all$pt_variant <- paste(vaf.all$patient, vaf.all$variant, sep='_')
sort(table(vaf.all[!duplicated(vaf.all$pt_variant), "variant"]), decreasing = T)
sort(table(vaf.all[!duplicated(vaf.all$pt_variant), "gene"]), decreasing = T)

# plot for a specific gene
variant <- 'JAK2_V617F'
vaf <- vaf.all[vaf.all$variant==variant, ]

# subset to patients that have multiple followup samples
vaf <- vaf[vaf$fo_tumor >0, ]
vaf$patient_tp <- paste(vaf$patient, vaf$tp_numeric, sep='_')
# take average at each time point
vaf.avg <- aggregate(vaf, by= list(vaf$patient_tp), FUN = function(x) {
    if(is.numeric(x)){
        mean(x)
    } else {
        x[1]
    }})

ggplot(vaf.avg, aes(x=tp_numeric, y=vaf, color=patient)) +
    geom_line() +
    theme_bw() +
    labs(title='Variant allele frequency: JAK2 V617F', subtitle='Mean value per timepoint, relative to diagnosis', x='Time (days) relative to diagnosis', y='Variant allele frequency')

# scale this relative to dx
vaf.scale <- vaf.avg
first.vaf <- vaf.avg[vaf.avg$tp_numeric==0, 'vaf']
names(first.vaf) <- vaf.avg[vaf.avg$tp_numeric==0, 'patient']
vaf.scale$vaf <- (vaf.scale$vaf - first.vaf[vaf.scale$patient]) / first.vaf[vaf.scale$patient]

p1 <- ggplot(vaf.scale, aes(x=tp_numeric, y=vaf, color=patient)) +
    geom_line(size=0.9) +
    theme_bw() +
    scale_y_continuous(labels=scales::percent, limits = c(-0.75, 0.25)) +
    geom_hline(yintercept = 0, lty=2, col='grey50') +
    labs(title='Variant allele frequency: JAK2 V617F', subtitle='Mean value per timepoint, relative to diagnosis', x='Time (days) relative to diagnosis', y='Variant allele frequency relative to diagnosis')
p1

# some sort of label for if patients are increasing or decreasing
major.decrease <- -0.25
minor.decrease <- 0
minor.increase <- 0.25
# value for last sample from a patient
vaf.final <- vaf.scale[order(vaf.scale$tp_numeric, decreasing = T), ]
vaf.final <- vaf.final[!duplicated(vaf.final$patient),]
vaf.final$classification <- 'minor decrease'
vaf.final[vaf.final$vaf < major.decrease, "classification"] <- 'major decrease'
vaf.final[vaf.final$vaf > 0 & vaf.final$vaf < minor.increase, "classification"] <- 'minor increase'
vaf.final[vaf.final$vaf > minor.increase, "classification"] <- 'major increase'
patient.classification <- vaf.final$classification
names(patient.classification) <- vaf.final$patient

vaf.scale$classification <- patient.classification[vaf.scale$patient]
vaf.export <- vaf.scale[,c("patient", "tp_numeric", "variant", 'vaf', 'classification')]
write.table(vaf.export, '~/pcloud_sync/project/vaf_decrease_plotting/JAK2_V617F_data_export.tsv', sep = '\t', quote=F, row.names = F, col.names = T)

pal = c('minor decrease'='lightblue', 'major decrease'='midnightblue', 'minor increase' = 'salmon', 'major increase' = 'darkred')
p2 <- ggplot(vaf.scale, aes(x=tp_numeric, y=vaf, group=patient,  color=classification)) +
    geom_line(size=0.9) +
    theme_bw() +
    scale_color_manual(values=pal) +
    scale_y_continuous(labels=scales::percent, limits = c(-0.75, 0.25)) +
    geom_hline(yintercept = 0, lty=2, col='grey50') +
    labs(title='Mutant allele frequency: JAK2 V617F', subtitle='Mean value per timepoint, relative to enrollment', x='Time (days) relative to enrollment', y='Mutant allele frequency relative to enrollment')
p2

pdf('~/pcloud_sync/project/vaf_decrease_plotting/JAK2_V617F_CTP_102_patients.pdf', height=5, width=6)
p2
p1
dev.off()
png('~/pcloud_sync/project/vaf_decrease_plotting/JAK2_V617F_CTP_102_patients.png', height = 750, width = 1000,res = 175)
p2
dev.off()


#############################################################################
# 2021-05-08
# plot all variants, or at least all genes
gene.table <- sort(table(vaf.all[!duplicated(vaf.all$pt_variant), "gene"]), decreasing = T)
plot.genes <- names(gene.table)[gene.table >1]
major.decrease <- -0.25
minor.decrease <- 0
minor.increase <- 0.25
pal = c('minor decrease'='lightblue', 'major decrease'='midnightblue', 'minor increase' = 'salmon', 'major increase' = 'darkred')



# pdf('~/pcloud_sync/project/vaf_decrease_plotting/all_gene_decrease_plot.pdf', height=5, width=6)
for (gene in plot.genes){
    print(gene)
    vaf <- vaf.all[vaf.all$gene==gene, ]
    vaf$patient_tp <- paste(vaf$patient, vaf$tp_numeric, sep='_')
    # take average at each time point
    vaf.avg <- aggregate(vaf, by= list(vaf$patient_tp), FUN = function(x) {
        if(is.numeric(x)){
            mean(x)
        } else {
            x[1]
        }})
    # scale this relative to dx
    vaf.scale <- vaf.avg
    first.vaf <- vaf.avg[vaf.avg$tp_numeric==0, 'vaf']
    names(first.vaf) <- vaf.avg[vaf.avg$tp_numeric==0, 'patient']
    vaf.scale$vaf <- (vaf.scale$vaf - first.vaf[vaf.scale$patient]) / first.vaf[vaf.scale$patient]

    vaf.final <- vaf.scale[order(vaf.scale$tp_numeric, decreasing = T), ]
    vaf.final <- vaf.final[!duplicated(vaf.final$patient),]
    vaf.final$classification <- 'minor decrease'
    vaf.final[vaf.final$vaf < major.decrease, "classification"] <- 'major decrease'
    vaf.final[vaf.final$vaf > 0 & vaf.final$vaf < minor.increase, "classification"] <- 'minor increase'
    vaf.final[vaf.final$vaf > minor.increase, "classification"] <- 'major increase'
    patient.classification <- vaf.final$classification
    names(patient.classification) <- vaf.final$patient
    vaf.scale$classification <- patient.classification[vaf.scale$patient]
    # So V617F appears on plot
    if (gene=='JAK2'){
        gene <- 'JAK2 V617F'
    }
    p2 <- ggplot(vaf.scale, aes(x=tp_numeric, y=vaf, group=patient,  color=classification)) +
        geom_line(size=0.9) +
        theme_bw() +
        scale_color_manual(values=pal) +
        scale_y_continuous(labels=scales::percent) +
        geom_hline(yintercept = 0, lty=2, col='grey50') +
        labs(title=paste('Mutant allele frequency:', gene),
             subtitle='Mean value per timepoint, relative to enrollment',
             x='Time (days) relative to enrollment',
             y='Mutant allele frequency relative to enrollment sample')
    # do some png exports instead
    outf <- paste0('~/pcloud_sync/project/vaf_decrease_plotting/png/maf_', gsub(' ', '_', gene), ".png")
    png(outf, height = 750, width = 1000,res = 175)
    print(p2)
    dev.off()
}

# dev.off()
