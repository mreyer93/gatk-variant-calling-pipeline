# combine annotations with a filtered VCF file to produce a final report
# this version specific to the TvN files

# if launched from snakemake set args here
if (exists('snakemake')) {
    vcf.input.f <- snakemake@input[['vcf_funcotator']]
    cadd.input.f <- snakemake@input[['vcf_cadd']]
    vcf.output.unfilt.f <- snakemake@output[['vcf_unfilt']]
    vcf.output.filt.f <- snakemake@output[['vcf_filt']]
    vaf.min <- snakemake@params[['filter_min_af']]
    depth.min <- snakemake@params[['filter_min_dp']]
    min.cadd <- snakemake@params[['filter_min_cadd']]
    keep.gatk.filters <- snakemake@params[['filter_GATK']]

} else{
    # parse args from command line
    args = commandArgs(trailingOnly=TRUE)
    if (length(args)!=4) {
        stop("Four arguments must be supplied: funcotator_input cadd_input unfiltered_output filtered_output", call.=FALSE)
    }
    vcf.input.f <- args[1]
    cadd.input.f <- args[2]
    vcf.output.unfilt.f <- args[3]
    vcf.output.filt.f <- args[4]
    vaf.min <- 0.02
    depth.min <- 5
    keep.gatk.filters <- c('PASS', 'normal_artifact')
}

# testing arguments
if(F){
    vcf.input.f <- '/path/to/data/pipeline_test/ctp_101_test_output/03_variants_TvN/02_variant_annotations/annotate_Funcotator/001-002_Dx_filtered_fucnc.vcf'
    cadd.input.f <- '/path/to/data/pipeline_test/ctp_101_test_output/03_variants_TvN/02_variant_annotations/annotate_CADD/simple_table/001-002_Dx_filtered_CADD_simple.vcf'
    vcf.output.unfilt.f <- '/path/to/data/pipeline_test/ctp_101_test_output/03_variants_TvN/02_variant_annotations/annotations_combined/unfiltered/001-002_Dx_unfiltered_annotated.vcf'
    vcf.output.filt.f <- '/path/to/data/pipeline_test/ctp_101_test_output/03_variants_TvN/02_variant_annotations/annotations_combined/filtered/001-002_Dx_filtered_annotated.vcf'
    vaf.min <- 0.02
    depth.min <- 5
    keep.gatk.filters <- c('PASS', 'normal_artifact')
}

# read once to get position of #CHROM
vcf <- read.table(vcf.input.f, sep='\t', quote='', comment.char = '', header = F, fill=T, nrows = 10000)
start.line <- which(vcf$V1=='#CHROM')
vcf <- read.table(vcf.input.f, sep='\t', skip = (start.line-1), quote='', comment.char = '', header = T, fill=T, check.names = F)
colnames(vcf)[1] <- 'CHROM'

# catch case with no variants
if(nrow(vcf)==0){
    write.table(vcf, vcf.output.unfilt.f, sep='\t', quote=F, row.names = F, col.names = T)
    write.table(vcf, vcf.output.filt.f, sep='\t', quote=F, row.names = F, col.names = T)

} else {
    # get the info from the sample columns we actually want to keep
    sample.cols <- 10:ncol(vcf)
    sample.df <- vcf[, sample.cols, drop=F]
    # apply across each of these to keep only the allelic fraction and depth
    names <- colnames(sample.df)
    nsamples <- length(names)
    sample.af.df <- as.data.frame(matrix(apply(sample.df,1:2, function(x) strsplit(x, split = ':')[[1]][c(3,4)]), ncol=nsamples*2, byrow = T))
    colnames(sample.af.df) <- paste0(rep(c('AF_', 'DP_'),times=nsamples), names)

    # apply variant allele frequency filter
    # must exceed this AF in at least one tumor sample
    # multiallelics give me trouble here
    # should check any of the multi are greater?
    af.columns <- (1:nsamples) *2 -1
    dp.columns <- (1:nsamples) *2
    keep.vars <- apply(sample.af.df[,af.columns, drop=F], 1, function(x) {
        any(sapply(x, function(y) { # check for multialleics
            if (grepl(',', y)){
                multi.afs <- strsplit(y, split=',')[[1]]
                return(any(as.numeric(multi.afs) >= vaf.min))
            } else {
                return(as.numeric(y) >= vaf.min)
            }}))
        })
    # calculate the depth filter
    keep.dp <- apply(sample.af.df[, dp.columns, drop=F], 1, function(x) any(x >= depth.min))

    # FUNCOTATION parsing
    # INFO string, we just want some info on the gene change
    funcotation.str <- sapply(vcf$INFO, function(x) {
        a <- strsplit(x, split=';', useBytes = TRUE)[[1]]
        to.ret <- a[which(sapply(a, function(b) substr(b, 1, 11)) == 'FUNCOTATION')]
        to.ret <- gsub('FUNCOTATION=\\[','',to.ret)
        to.ret <- gsub('\\]','',to.ret)
        to.ret <- gsub('\\"','',to.ret)
        return(to.ret)
    })
    names(funcotation.str) <- NULL
    keep.func <- c(1,2,6,7,8,19, 12,17,18)
    funcotation.list <- lapply(funcotation.str, function(x) strsplit(x, split='|', fixed=T, useBytes = TRUE)[[1]][keep.func])
    funcotation.df.keep <- do.call(rbind, funcotation.list)
    rownames(funcotation.df.keep) <- NULL
    new.colnames <- c('Gencode_34_hugoSymbol', 'Gencode_34_ncbiBuild', 'Gencode_34_variantClassification', 'Gencode_34_secondaryVariantClassification',
                   'Gencode_34_variantType', 'Gencode_34_proteinChange', 'Gencode_34_genomeChange', 'Gencode_34_cDnaChange','Gencode_34_codonChange')
    colnames(funcotation.df.keep) <- new.colnames
    colnames(funcotation.df.keep) <- gsub('Gencode_34_', '', colnames(funcotation.df.keep))

    # add this back to the main df
    vcf.func <- cbind(vcf[, c('CHROM', 'POS', 'REF', 'ALT', 'FILTER')], funcotation.df.keep, sample.af.df)
    rm(vcf)

    # get CADD score information
    cadd <- read.table(cadd.input.f, sep='\t', quote='', comment.char = '', header = T, fill=T, skip = 1)
    colnames(cadd)[1] <- 'CHROM'
    # unique names for cadd variants and other variants
    cadd$uniq.var <- paste(gsub('chr', '', cadd$CHROM), cadd$Pos, cadd$Ref, cadd$Alt, sep='_')
    # remove duplicates here because all I care about is the phred score
    cadd <- cadd[!duplicated(cadd$uniq.var), ]
    rownames(cadd) <- cadd$uniq.var
    vcf.uniq.var <- paste(gsub('chr', '', vcf.func$CHROM), vcf.func$POS, vcf.func$REF, vcf.func$ALT, sep='_')
    vcf.func$CADD_phred <- cadd[vcf.uniq.var, "PHRED"]
    keep.cadd <- vcf.func$CADD_phred >= min.cadd
    # change col names so cadd is earlier
    col.order <- c(1:14, ncol(vcf.func), 15:(ncol(vcf.func)-1))
    vcf.func <- vcf.func[, col.order]
    
    # apply the VAF filter and depth filter
    vcf.func.filt <- vcf.func[keep.vars & keep.dp & keep.cadd, ]
    # apply GATK filter
    vcf.func.filt <- vcf.func.filt[vcf.func.filt$FILTER %in% keep.gatk.filters, ]
    # sample.af.df.filt <- sample.af.df[keep.vars & keep.dp, ]

    # sort by decreasing CADD_phred
    vcf.func <- vcf.func[order(vcf.func$CADD_phred, decreasing=T) ,]
    vcf.func.filt <- vcf.func.filt[order(vcf.func.filt$CADD_phred, decreasing=T) ,]

    # write out one version with the filter applied, and one without
    write.table(vcf.func, vcf.output.unfilt.f, sep='\t', quote=F, row.names = F, col.names = T)
    write.table(vcf.func.filt, vcf.output.filt.f, sep='\t', quote=F, row.names = F, col.names = T)
}
