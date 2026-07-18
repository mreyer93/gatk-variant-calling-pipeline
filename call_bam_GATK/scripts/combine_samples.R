# combine annotated variant calls from multiple samples from the same patient
# by looping through each patient in the metadata file
# this is going to assume that decomposition is not necessary, as all the variants
# we expect among this set should be simple
# JUST FOR THE TUMOR VS REFERENCE CALLS

# input files from snakemake
metadata.f <- snakemake@input[['metadata']]
annot_combined_dir <- snakemake@params[['annot_combined_dir']]
outdir.base <- snakemake@params[['outdir_base']]
# filtering args
filter.min.cadd <- as.numeric(snakemake@params[['filter_min_cadd']])
filter.min.af <- as.numeric(snakemake@params[['filter_min_af']])
filter.min.dp <- as.numeric(snakemake@params[['filter_min_dp']])
keep.gtak.filters <- snakemake@params[['filter_GATK']]

# TESTING ARGS
# temp processing outside of snakemake
if(F){
    # EGA data stuff
    metadata.f <- '~/local_data/EGA_data/bam_metadata_uniq_samp.tsv'
    annot_combined_dir <- '~/local_data/EGA_data/02_variants_reference/annotations_combined/'
    outdir.base <- '~/local_data/EGA_data/02_variants_reference/new_scores/'
    filter.min.cadd <- 20
    filter.min.af <- 0.02
    filter.min.dp <- 5
    # EGA data EXOME
    metadata.f <- '~/local_data/EGA_data/bam_metadata_uniq_samp.tsv'
    annot_combined_dir <- '~/local_data/EGA_data/02_variants_reference_exome/annotations_combined/'
    outdir.base <- '~/local_data/EGA_data/02_variants_reference_exome/new_scores/'
    filter.min.cadd <- 20
    filter.min.af <- 0.02
    filter.min.dp <- 5
    # CTP-102
    metadata.f <- '~/local_data/data_store/ctp_102_calling/metadata.tsv'
    annot_combined_dir <- '~/local_data/data_store/ctp_102_calling/output/02_variants_reference/annotations_combined/'
    outdir.base <- '~/local_data/data_store/ctp_102_calling/output/02_variants_reference/new_scores/'
    filter.min.cadd <- 20
    filter.min.af <- 0.02
    filter.min.dp <- 5
    # CTP-101
    metadata.f <- '~/local_data/data_store/ctp_101_calling/bam_metadata_all.tsv'
    annot_combined_dir <- '~/local_data/data_store/ctp_101_calling/output/02_variants_reference/annotations_combined/'
    outdir.base <- '~/local_data/data_store/ctp_101_calling/output/02_variants_reference/new_scores/'
    filter.min.cadd <- 20
    filter.min.af <- 0.02
    filter.min.dp <- 5
}
# input directories
indir.unfilt <- file.path(annot_combined_dir, 'unfiltered')
# output directories
outdir.unfilt <- file.path(outdir.base, 'unfiltered')
outdir.filt <- file.path(outdir.base, 'filtered')
outdir.filt.gatk <- file.path(outdir.base, 'filtered_GATK')
outdir.bedfiles <- file.path(outdir.base, 'to_compute_bedfiles')
dir.create(outdir.unfilt, recursive = T, showWarnings = F)
dir.create(outdir.filt, recursive = T, showWarnings = F)
dir.create(outdir.filt.gatk, recursive = T, showWarnings = F)
dir.create(outdir.bedfiles, recursive = T, showWarnings = F)
# read metadata
metadata <- read.table(metadata.f, sep='\t', header=T, comment.char ='', fill=T, quote='')
rownames(metadata) <- metadata$sample
pts <- unique(metadata$patient)

# function to take mean value of allelic fractions for mutliallelic sites
# operates on a whole data.fame
multi_to_mean_mat <- function(df){
    new.df <- apply(df, 1:2, function(x) {
        if(length(grep(',', x)>0)){
            return(mean(as.numeric(strsplit(x, ",")[[1]], na.rm = T)))
        } else{
            return(as.numeric(x))
        }
    })
    return(as.matrix(new.df))
}

# loop over each patient
for (pt in pts){
    print(pt)
    samples <- sort(metadata[metadata$patient==pt, "sample"])
    sample.colnames <- c(sapply(samples, function(x) paste0(c('AF_', 'DP_'), x)))
    tumor.samples <- samples[metadata[samples, "tumor_normal"] =='T']
    normal.samples <- samples[metadata[samples, "tumor_normal"] =='N']
    tumor.sample.colnames <- c(sapply(tumor.samples, function(x) paste0(c('AF_', 'DP_'), x)))
    normal.sample.colnames <- c(sapply(normal.samples, function(x) paste0(c('AF_', 'DP_'), x)))

    # get vcf files assuming they follow my convention
    vcf.unfilt.f <- sapply(samples, function(x) file.path(indir.unfilt, paste0(x, '_unfiltered_annotated.vcf')))
    if(!all(file.exists(vcf.unfilt.f))){
        print(vcf.unfilt.f)
        stop('Some vcf files do not exist')
    }

    # read in the dataframes of vcf files
    outf.unfilt <- file.path(outdir.unfilt, paste0(pt, '.vcf'))
    outf.filt <- file.path(outdir.filt, paste0(pt, '.vcf'))
    outf.filt.gatk <- file.path(outdir.filt.gatk, paste0(pt, '.vcf'))
    df.list <- lapply(vcf.unfilt.f, function(x) {
        a <- read.table(x, header=T, sep='\t', quote='', comment.char = '', check.names = F)
        a <- a[order(a$CHROM, a$POS, a$REF, a$ALT),]
        a$key <- paste(a$CHROM, a$POS, a$REF, a$ALT, sep='_')
        return(a)
    })

    df.annot <- do.call(rbind, lapply(df.list, function(x) x[,colnames(x)[!colnames(x) %in% sample.colnames]]))
    df.annot <- df.annot[!duplicated(df.annot$key), ]
    rownames(df.annot) <- df.annot$key

    # take only key and AF, DP
    df.list.short <- lapply(samples, function(x) df.list[[x]][,c('key',  paste0(c('AF_', 'DP_'), x))])

    # merge them all on the same key...
    df.all <- Reduce(function(...) merge(..., by='key', all=TRUE, sort=F, no.dups=T), df.list.short)
    df.all <- df.all[!duplicated(df.all$key), ]
    rownames(df.all) <- df.all$key

    # combine with annotations
    df.comb <- cbind(df.annot, df.all[rownames(df.annot), sample.colnames, drop=F])
    # sort by decreasing CADD_phred
    df.comb <- df.comb[order(df.comb$CADD_phred, decreasing = T),]
    # add gene_protein annotation
    df.comb$proteinChange[is.na(df.comb$proteinChange)] <- ''
    df.comb$gene_protein <- paste(df.comb$hugoSymbol, df.comb$proteinChange, sep='_')
    df.comb$gene_protein[df.comb$proteinChange==''] <- paste(df.comb$hugoSymbol[df.comb$proteinChange==''],
                                                             df.comb$variantClassification[df.comb$proteinChange==''], sep='_')

    # change columns to be tumor then normal
    annot.cols <- c(colnames(df.comb)[1:15], 'gene_protein')
    df.comb <- df.comb[, unlist(c(annot.cols, tumor.sample.colnames, normal.sample.colnames))]

    # some metric of tumor vs normal differece
    # could do a stat test, anova, t test, something here
    # limited sample number is going to be the problem
    # tumor vs normal mean log2 fold change?
    # two versions, one where we count NA as zero, and one where it's ignored
    # can only do this if we have at least one tumor sample
    if (length(tumor.samples)>0){
        T.mat <- multi_to_mean_mat(df.comb[, paste0('AF_', tumor.samples), drop=F])
        T.mat.na0 <- T.mat
        T.mat.na0[is.na(T.mat.na0)] <- 0
        mean.T <- rowMeans(T.mat.na0)
        mean.T.nna <- rowMeans(T.mat, na.rm = T)
        # make a version for depth too
        T.mat.dp <- multi_to_mean_mat(df.comb[, paste0('DP_', tumor.samples), drop=F])

        # add mean to the df
        df.comb$mean.T.nna <- round(mean.T.nna, 3)
    } else{
        df.comb$mean.T.nna <- NA
    }

    # can only do this if we have at least one normal sample
    if (length(normal.samples)>0){
        N.mat <- multi_to_mean_mat(df.comb[, paste0('AF_', normal.samples), drop=F])
        N.mat.na0 <- N.mat
        N.mat.na0[is.na(N.mat.na0)] <- 0
        mean.N <- rowMeans(N.mat.na0)
        mean.N.nna <- rowMeans(N.mat, na.rm = T)
        # add mean to the df
        df.comb$mean.N.nna <- round(mean.N.nna, 3)
    } else {
        df.comb$mean.N.nna <- NA
    }

    # can only do this if we have both tumor and normal samples
    if ((length(normal.samples)>0) & (length(tumor.samples)>0)){
        # log fold change calcs
        # epsilon for division
        e <- 0.001
        # cap at 2.5
        max.lfc <- 2.5
        l2fc <- log2(mean.T/(mean.N + e))
        l2fc[l2fc > max.lfc] <- max.lfc
        l2fc.nna <- log2(mean.T.nna/(mean.N.nna + e))
        l2fc.nna[l2fc.nna > max.lfc] <- max.lfc
        # this metric seems good, which combines CADD score, log2fc and mean tumor AF
        caddlfc <- df.comb$CADD_phred * l2fc * mean.T
        caddlfc.nna <- df.comb$CADD_phred * l2fc.nna * mean.T.nna

        decide.score <- apply(data.frame(round(caddlfc, 2), round(caddlfc.nna, 2)), 1, function(x) max(x, na.rm = T))
        df.comb$score.new <- decide.score
    } else {
        df.comb$score.new <- df.comb$CADD_phred
    }
    # sort by new score
    df.comb <- df.comb[order(df.comb$score.new, decreasing = T),]

    # change column orders to have samples last
    new.colorder <- unlist(c(annot.cols, 'mean.T.nna', 'mean.N.nna', 'score.new', tumor.sample.colnames, normal.sample.colnames))
    df.comb <- df.comb[, new.colorder]
    # write out unfiltered
    write.table(df.comb, outf.unfilt, sep='\t', quote=F, row.names = F, col.names = T)

    # filtered version
    # just cadd score and AF filters on tumor samples
    # look at T.mat to see if any are above the threshold
    keep.af <- apply(T.mat, 1, function(x) any(x >= filter.min.af, na.rm = T))[rownames(df.comb)]
    # look at depth in tumor samples
    keep.dp <- apply(T.mat.dp, 1, function(x) any(x >= filter.min.dp, na.rm = T))[rownames(df.comb)]

    df.comb.filt <- df.comb[!is.na(df.comb$CADD_phred) & df.comb$CADD_phred >= filter.min.cadd & keep.af & keep.dp, ]
    # write out filtered
    write.table(df.comb.filt, outf.filt, sep='\t', quote=F, row.names = F, col.names = T)
    # and another version passing al GATK filters``
    write.table(df.comb.filt[df.comb.filt$FILTER %in% keep.gtak.filters,], outf.filt.gatk, sep='\t', quote=F, row.names = F, col.names = T)

    # get the locations that are NA that we want to genotype in the other samples
    # or at least get the depth for
    # one for each sample?
    # and write out a bed file for this sample
    # this really only applies to normal samples because we selected variants with
    # a min depth in the tumor
    if (length(normal.samples)>0 & length(tumor.samples) > 0){
        for (s in normal.samples){
            out.bed <- file.path(outdir.bedfiles, paste0(s,'.bed'))
            this.col <- paste0('DP_', s)
            check.df <- df.comb.filt[is.na(df.comb.filt[,this.col]),]
            check.df <- check.df[!is.na(check.df$CHROM),]
            check.df$key <- paste(check.df$CHROM, check.df$POS, sep='_')
            check.df <- check.df[order(check.df$key),]
            bed.df <- data.frame(chr=check.df$CHROM, start=check.df$POS-1, stop=check.df$POS)
            write.table(bed.df, out.bed, sep='\t', row.names = F, col.names = F, quote=F)
        }
    }
}
