# how well are exons covered in BOD1L1 in EGA data?
library(ggplot2)
library(reshape2)
options(stringsAsFactors = F)
cov <- read.table('~/data/EGA_data/01_prepare_bam/BOD1L1_coverage.txt', sep='\t', quote='', header=F)
cov.df <- data.frame(exon=cov$V5, coverage = cov$V10)
cov.df$exon <- sapply(cov.df$exon, function(x) strsplit(x, split='_')[[1]][3])
cov.df <- cov.df[order(as.numeric(cov.df$exon)), ]
cov.df$exon <- factor(cov.df$exon, levels <- unique(cov.df$exon))

# truncate this to 600
max.cov <- 600
cov.df$coverage[cov.df$coverage > max.cov] <- max.cov

ggplot(cov.df, aes(y=coverage, x=exon)) +
    geom_boxplot() +
    theme_bw() +
    labs(title='Coverage of BOD1L1 in EGA')

# do the same thing for CTP 102 as comparison
cov <- read.table('~/data/ctp_102_calling/gatk_calling_bam/output/01_prepare_bam/BOD1L1_coverage.tsv', sep='\t', quote='', header=F)
cov.df <- data.frame(exon=cov$V5, coverage = cov$V10)
cov.df$exon <- sapply(cov.df$exon, function(x) strsplit(x, split='_')[[1]][3])
cov.df <- cov.df[order(as.numeric(cov.df$exon)), ]
cov.df$exon <- factor(cov.df$exon, levels <- unique(cov.df$exon))

# truncate this to 600
max.cov <- 2000
cov.df$coverage[cov.df$coverage > max.cov] <- max.cov

ggplot(cov.df, aes(y=coverage, x=exon)) +
    geom_boxplot() +
    theme_bw() +
    labs(title='Coverage of BOD1L1 in CTP_102')


########################################################################################################
# do the same thing but for MAP1B ######################################################################
########################################################################################################
cov <- read.table('~/data/EGA_data/01_prepare_bam/MAP1B_coverage.txt', sep='\t', quote='', header=F)
cov.df <- data.frame(exon=cov$V5, coverage = cov$V10)
cov.df$exon <- sapply(cov.df$exon, function(x) strsplit(x, split='_')[[1]][3])
cov.df <- cov.df[order(as.numeric(cov.df$exon)), ]
cov.df$exon <- factor(cov.df$exon, levels <- unique(cov.df$exon))

# truncate this to 600
max.cov <- 600
cov.df$coverage[cov.df$coverage > max.cov] <- max.cov

ggplot(cov.df, aes(y=coverage, x=exon)) +
    geom_boxplot() +
    theme_bw() +
    labs(title='Coverage of MAP1B in EGA')

# do the same thing for CTP 102 as comparison
cov <- read.table('~/data/ctp_102_calling/gatk_calling_bam/output/01_prepare_bam/MAP1B_coverage.tsv', sep='\t', quote='', header=F)
cov.df <- data.frame(exon=cov$V5, coverage = cov$V10)
cov.df$exon <- sapply(cov.df$exon, function(x) strsplit(x, split='_')[[1]][3])
cov.df <- cov.df[order(as.numeric(cov.df$exon)), ]
cov.df$exon <- factor(cov.df$exon, levels <- unique(cov.df$exon))

# truncate this to 600
max.cov <- 2000
cov.df$coverage[cov.df$coverage > max.cov] <- max.cov

ggplot(cov.df, aes(y=coverage, x=exon)) +
    geom_boxplot() +
    theme_bw() +
    labs(title='Coverage of MAP1B in CTP_102')
