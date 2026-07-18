# This script re-computes the new scores with the new genotypes
# after they've been re-computed for the normal samples
# at the specific missing sites
metadata.f <- snakemake@input[['metadata']]
indir.base <- snakemake@params[['indir_base']]
indir.vcf <- snakemake@params[['indir_vcf']]
outdir.base <- snakemake@params[['outdir_base']]
# filtering args
filter.min.cadd <- as.numeric(snakemake@params[['filter_min_cadd']])
filter.min.af <- as.numeric(snakemake@params[['filter_min_af']])
filter.min.dp <- as.numeric(snakemake@params[['filter_min_dp']])
keep.gtak.filters <- snakemake@params[['filter_GATK']]

# TESTING ARGS
if(F){
    # CTP DATA
    trial <- 'ctp_102'
    metadata.f <- paste0('~/local_data/data_store/', trial, '_calling/bam_metadata.tsv')
    indir.base <- paste0('~/local_data/data_store/', trial, '_calling/output/02_variants_reference/new_scores/filtered')
    indir.vcf <- paste0('~/local_data/data_store/', trial, '_calling/output/02_variants_reference/new_scores/computed_vcfs')
    outdir.base <- paste0('~/local_data/data_store/', trial, '_calling/output/02_variants_reference/REDO_new_scores/')
    filter.min.cadd <- 20
    filter.min.af <- 0.02
    filter.min.dp <- 5
    # EGA data EXOME
    metadata.f <- '~/local_data/EGA_data/bam_metadata_uniq_samp.tsv'
    indir.base <- '~/local_data/EGA_data/02_variants_reference_exome/new_scores/filtered'
    indir.vcf <- '~/local_data/EGA_data/02_variants_reference_exome/new_scores/computed_vcfs'
    outdir.base <- '~/local_data/EGA_data/02_variants_reference_exome/REDO_new_scores/'
    filter.min.cadd <- 20
    filter.min.af <- 0.02
    filter.min.dp <- 5
}
# output directories
outdir.filt <- file.path(outdir.base, 'filtered')
outdir.filt.gatk <- file.path(outdir.base, 'filtered_GATK')
dir.create(outdir.filt, recursive = T, showWarnings = F)
dir.create(outdir.filt.gatk, recursive = T, showWarnings = F)
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
    # read in the df.comb.filt for this patient
    inf.filt <- file.path(indir.base, paste0(pt, '.vcf'))
    df.comb.filt <- read.table(inf.filt, header=T, sep='\t', quote='',check.names = F)
    # output files
    outf.filt <- file.path(outdir.filt, paste0(pt, '.vcf'))
    outf.filt.gatk <- file.path(outdir.filt.gatk, paste0(pt, '.vcf'))

    # can only do this if we have both tumor and normal samples
    if ((length(normal.samples)>0) & (length(tumor.samples)>0)){
        for (s in normal.samples){
            s.colnames <- paste0(c('AF_', 'DP_'), s)
            df.comb.filt$key <- paste(df.comb.filt$CHROM, df.comb.filt$POS, sep='_')
            dup.keys <- df.comb.filt$key[duplicated(df.comb.filt$key)]
            if (length(dup.keys)>0){
                # print(df.comb.filt[df.comb.filt$key %in% dup.keys, ])
                dup.keys <- paste(dup.keys, 1:length(dup.keys), sep='_')
                df.comb.filt$key[duplicated(df.comb.filt$key)] <- dup.keys
            }
            rownames(df.comb.filt) <- df.comb.filt$key
            vcf.input.f <- file.path(indir.vcf, paste0(s,'.vcf'))
            vcf <- read.table(vcf.input.f, sep='\t', quote='', comment.char = '', header = F, fill=T, nrows = 10000)
            start.line <- which(vcf$V1=='#CHROM')
            vcf <- read.table(vcf.input.f, sep='\t', skip = (start.line-1), quote='', comment.char = '', header = T, fill=T, check.names = F)
            colnames(vcf)[1] <- 'CHROM'
                # only if theres data in the vcf file
                if(nrow(vcf) >0){
                    vcf$key <- paste(vcf$CHROM, vcf$POS, sep='_')
                    vcf <- vcf[order(vcf$key),]
                    # depth info is in the first info field
                    vcf$depth <- sapply(vcf$INFO, function(x) as.integer(gsub('DP=', '', grep('DP=', strsplit(x, split=';')[[1]],value=T))))
                    # only first alt considered here
                    vcf$AF <- sapply(vcf$INFO, function(x) as.numeric(gsub('AF1=', '', grep('AF1=', strsplit(x, split=';')[[1]],value=T)[1])))
                    # remove duplicate entries without an alternate allele
                    vcf <- vcf[order(vcf$ALT, decreasing = T),]
                    vcf <- vcf[!(duplicated(vcf$key)),]
                    rownames(vcf) <- vcf$key

                    # set all the ones with NO ALT
                    vcf.noalt <- vcf[vcf$ALT=='.', ]
                    df.comb.filt[vcf.noalt$key, s.colnames] <- vcf.noalt[,c("AF", "depth")]

                    vcf.alt <- vcf[vcf$ALT!='.', ]
                    for (key in vcf.alt$key){
                        tumor.alt <- df.comb.filt[key, "ALT"]
                        this.alt <- vcf.alt[key, "ALT"]
                        if(tumor.alt == this.alt){
                            df.comb.filt[key, s.colnames] <- vcf.alt[key,c("AF", "depth")]
                        } else {
                            df.comb.filt[key, s.colnames] <- c(NA, vcf.alt[key,c("depth")])
                        }

                }
            } else {
                df.comb.filt[key, s.colnames] <- NA
            }

            # so only a few na regions should remain
            # still.na <- df.comb.filt[is.na(df.comb.filt[,s.colnames[2]]), ]
        }


        # after all normal samples have their data added, we can take the set of non-na keys
        # df.comb.filt.fixed <- df.comb.filt[!is.na(df.comb.filt[, s.colnames[2]])]

        # re-compute the N.mean and the new.scores
        N.mat <- multi_to_mean_mat(df.comb.filt[, paste0('AF_', normal.samples), drop=F])
        N.mat.na0 <- N.mat
        N.mat.na0[is.na(N.mat.na0)] <- 0
        mean.N <- rowMeans(N.mat.na0)
        mean.N.nna <- rowMeans(N.mat, na.rm = T)
        # add mean to the df
        df.comb.filt$mean.N.nna <- round(mean.N.nna, 3)

        # can only do this if we have both tumor and normal samples
        if ((length(normal.samples)>0) & (length(tumor.samples)>0)){
            # log fold change calcs
            # epsilon for division
            e <- 0.001
            # cap at 2.5
            max.lfc <- 2.5
            # l2fc <- log2(mean.T/(mean.N + e))
            # l2fc[l2fc > max.lfc] <- max.lfc
            l2fc.nna <- log2(df.comb.filt$mean.T.nna/(mean.N.nna + e))
            l2fc.nna[l2fc.nna > max.lfc] <- max.lfc
            # this metric seems good, which combines CADD score, log2fc and mean tumor AF
            # caddlfc <- df.comb$CADD_phred * l2fc * mean.T
            caddlfc.nna <- df.comb.filt$CADD_phred * l2fc.nna * df.comb.filt$mean.T.nna

            df.comb.filt$score.new <- caddlfc.nna
        }

        # to pass here, need at least depth filter.min.dp in at least one normal sample
        normal.sample.dpnames <- c(sapply(normal.samples, function(x) paste0(c('DP_'), x)))
        keep.normal.depth <- apply(df.comb.filt[, normal.sample.dpnames, drop=F], 1, function(x) {
            if(all(is.na(x))){
                return(F)
            } else {
                return(max(x, na.rm=T) >= filter.min.dp)
            }
        })

        df.comb.filt <- df.comb.filt[keep.normal.depth, ]
        # remove key column
        df.comb.filt <- df.comb.filt[, !(colnames(df.comb.filt) %in% c('key'))]
        # sort by new score
        df.comb.filt <- df.comb.filt[order(df.comb.filt$score.new, decreasing = T),]
        df.comb.filt$score.new <- round(df.comb.filt$score.new, 2)
        write.table(df.comb.filt, outf.filt, sep='\t', quote=F, row.names = F, col.names = T)
        # and another version passing al GATK filters``
        write.table(df.comb.filt[df.comb.filt$FILTER %in% keep.gtak.filters,], outf.filt.gatk, sep='\t', quote=F, row.names = F, col.names = T)

    } else {
        # just write out the same file again????
        write.table(df.comb.filt, outf.filt, sep='\t', quote=F, row.names = F, col.names = T)
        # and another version passing al GATK filters``
        write.table(df.comb.filt[df.comb.filt$FILTER %in% keep.gtak.filters,], outf.filt.gatk, sep='\t', quote=F, row.names = F, col.names = T)

    }
}
