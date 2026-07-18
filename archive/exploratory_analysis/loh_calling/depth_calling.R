#!/usr/bin/env Rscript
# new LOH calling script on the project data
# working with the somatic variant calls, and depth from a single sample
library(ggplot2)
library(reshape2)
library(dplyr)

# working with the smoothed, joined depth file
depth_f <- '~/data/CTP_201/01_prepare_bam/depth_smoothed_join.tsv'
depth_df <- read.table(depth_f, sep='\t', quote='', header=T, check.names = F)
rownames(depth_df) <- depth_df$name
chr1_names <- depth_df[depth_df$chr=='chr1', "name"]
chr17_names <- depth_df[depth_df$chr=='chr17', "name"]
chr22_names <- depth_df[depth_df$chr=='chr22', "name"]

sample_names <- colnames(depth_df)[9:ncol(depth_df)]
depth_mat <- depth_df[,sample_names]

depth_mat_sample_norm <- apply(depth_mat, 2, function(x) x/mean(x,na.rm=T))

c1m <- melt(depth_mat_sample_norm[chr1_names,])

ggplot(c1m, aes(x=Var1, y=value, group=Var2)) + 
    geom_line()

# then norm the rows
depth_mat_both_norm <- t(apply(depth_mat_sample_norm, 1, function(x) x/mean(x,na.rm=T)))
c1m2 <- melt(depth_mat_both_norm[chr1_names,])
ggplot(c1m2, aes(x=Var1, y=value, group=Var2)) + 
    geom_line()
c17m2 <- melt(depth_mat_both_norm[chr17_names,])
c17m2$patient <- sapply(as.character(c17m2$Var2), function(x) strsplit(x, split="_")[[1]][1])
c17m2$patient_plot <- 'other'
c17m2$patient_plot[c17m2$patient == '055-204'] <- '055-204'
c17m2$patient_plot[c17m2$patient == '001-205'] <- '001-205'
c17m2$start <- depth_df[c17m2$Var1, "start"]
ggplot(c17m2, aes(x=start, y=log2(value), color=patient_plot)) + 
    geom_point(alpha=0.5)
ggplot(c17m2[c17m2$patient_plot %in% c('055-204','001-205'),], aes(x=Var1, y=value, group=Var2, color=Var2)) + 
    geom_line()
ggplot(c17m2[c17m2$patient_plot %in% c('055-204','001-205'),], aes(x=start, y=log2(value), group=Var2, color=Var2)) + 
    geom_point()

c22m2 <- melt(depth_mat_both_norm[chr22_names,])
c22m2$patient <- sapply(as.character(c22m2$Var2), function(x) strsplit(x, split="_")[[1]][1])
c22m2$patient_plot <- 'other'
c22m2$patient_plot[c22m2$patient == '034-206'] <- '034-206'
ggplot(c22m2, aes(x=Var1, y=value, group=Var2, color=patient_plot)) + 
    geom_point()


    # my assumptions
# seq depth is determined on a per-sample basis
    # so everything in a sample needs to get normed to some per-sample value
# coverage of a region is determined by seq-specific stuff
    # so that region will be consistently up or down depending on how things go


# in aggergate, this doesn't work too well, at least not on CTP-201 data. I can't 
# a priori predict the regions that are losses becuase there's so much noise.  