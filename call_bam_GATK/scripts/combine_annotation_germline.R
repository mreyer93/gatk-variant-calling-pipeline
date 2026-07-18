# combine annotations with a filtered VCF file to produce a final report

# if launched from snakemake set args here
# cadd.input.f and am.input.f are each optional - CADD and AlphaMissense annotation are
# both toggled independently via config (skip_cadd, alphamissense_file). AlphaMissense
# only scores missense SNVs - indels/nonsense/splice/non-coding variants will have NA
# in the AM_pathogenicity/AM_class columns even when it's enabled.
if (exists('snakemake')) {
    vcf.input.f <- snakemake@input[['vcf_funcotator']]
    cadd.input.f <- if (length(snakemake@input[['vcf_cadd']]) > 0) snakemake@input[['vcf_cadd']] else NULL
    am.input.f <- if (length(snakemake@input[['vcf_alphamissense']]) > 0) snakemake@input[['vcf_alphamissense']] else NULL
    vcf.output.f <- snakemake@output[['vcf']]
    depth.min <- snakemake@params[['filter_min_dp']]

} else{
    # parse args from command line
    args = commandArgs(trailingOnly=TRUE)
    if (length(args)!=4) {
        stop("Four arguments must be supplied: funcotator_input cadd_input unfiltered_output filtered_output", call.=FALSE)
    }
    vcf.input.f <- args[1]
    cadd.input.f <- args[2]
    am.input.f <- NULL
    vcf.output.f <- args[3]
    depth.min <- 5
}

# testing arguments
if(F){
    # multi-sample test with germline data from CTP_201
    vcf.input.f <- '~/data/CTP_201_SF3B1/07_joint_vcf/02_variant_annotations/funcotator.vcf'
    cadd.input.f <- '~/data/CTP_201_SF3B1/07_joint_vcf/02_variant_annotations/CADD_simple_table.vcf'
    vcf.output.f <- '~/data/CTP_201_SF3B1/07_joint_vcf/02_variant_annotations/annotations_combined.vcf'
    depth.min <- 5
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
    sample.names <- colnames(sample.df)
    
    # GT:AD:DP:GQ:PL
    data.mat <- apply(sample.df, 1:2, function(x) strsplit(x, split = ':')[[1]][1:5])
    #5D matrix with second and third axis corresponding to rows and columns (variants and samples)
    # and the first axis corresponding to GT, AD, etc
    dimnames(data.mat)[[1]] <- c('GT', 'AD', 'DP', 'GQ', 'PL')
    
    # mask uncalled variants with NA
    mask.na.gt <- data.mat[1,,]=="./." | data.mat[1,,]=="."
    data.mat[1,,][mask.na.gt] <- NA
    data.mat[2,,][mask.na.gt] <- NA
    data.mat[3,,][mask.na.gt] <- NA
    data.mat[4,,][mask.na.gt] <- NA
    data.mat[5,,][mask.na.gt] <- NA
    
    # apply filters based on depth
    # min depth default is 5, mask everything without that
    mask.na.depth <- data.mat[3,,] < depth.min
    data.mat[1,,][mask.na.depth] <- NA
    data.mat[2,,][mask.na.depth] <- NA
    data.mat[3,,][mask.na.depth] <- NA
    data.mat[4,,][mask.na.depth] <- NA
    data.mat[5,,][mask.na.depth] <- NA
    
    gt.mat <- data.mat[1,,]
    ad.mat <- data.mat[2,,]
    dp.mat <- apply(data.mat[3,,],2,as.numeric)
    gq.mat <- apply(data.mat[4,,],2,as.numeric)
    pl.mat <- data.mat[5,,]
    
    gt.mat[1:5,1:5]
    ad.mat[1:5,1:5]
    dp.mat[1:5,1:5]
    gq.mat[1:5,1:5]
    pl.mat[1:5,1:5]
    
    # convert GT to simple, biallelic
    sort(table(gt.mat), decreasing = T)
    gt.mat.simple <- gt.mat
    gt.mat.simple[gt.mat.simple == '0/0'] <- 0
    gt.mat.simple[gt.mat.simple == '0/1'] <- 1
    gt.mat.simple[gt.mat.simple == '0/2'] <- 0
    gt.mat.simple[gt.mat.simple == '1/2'] <- 1
    gt.mat.simple[gt.mat.simple == '1/1'] <- 2
    gt.mat.simple[gt.mat.simple == '0|0'] <- 0
    gt.mat.simple[gt.mat.simple == '0|1'] <- 1
    gt.mat.simple[gt.mat.simple == '1|1'] <- 2
    gt.mat.simple[gt.mat.simple == '1|2'] <- 1
    gt.mat.simple <- apply(gt.mat.simple,2,as.numeric)
    
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
    vcf.func <- cbind(vcf[, c('CHROM', 'POS', 'REF', 'ALT', 'FILTER')], funcotation.df.keep)
    # rm(vcf)
    
    vcf.uniq.var <- paste(gsub('chr', '', vcf.func$CHROM), vcf.func$POS, vcf.func$REF, vcf.func$ALT, sep='_')

    # get CADD score information, if CADD annotation was run
    if (!is.null(cadd.input.f)) {
        cadd <- read.table(cadd.input.f, sep='\t', quote='', comment.char = '', header = T, fill=T, skip = 1)
        colnames(cadd)[1] <- 'CHROM'
        # unique names for cadd variants and other variants
        cadd$uniq.var <- paste(gsub('chr', '', cadd$CHROM), cadd$Pos, cadd$Ref, cadd$Alt, sep='_')
        # remove duplicates here because all I care about is the phred score
        cadd <- cadd[!duplicated(cadd$uniq.var), ]
        rownames(cadd) <- cadd$uniq.var
        vcf.func$CADD_phred <- cadd[vcf.uniq.var, "PHRED"]
    }

    # get AlphaMissense score information, if AlphaMissense annotation was run.
    # Note: AlphaMissense only scores missense SNVs, so indels/nonsense/splice/non-coding
    # variants will be NA here even when this ran.
    if (!is.null(am.input.f)) {
        am <- read.table(am.input.f, sep='\t', quote='', comment.char = '', header = T, fill=T)
        # multiple transcripts can give multiple predictions per site - keep the first
        am$uniq.var <- paste(gsub('chr', '', am$CHROM), am$POS, am$REF, am$ALT, sep='_')
        am <- am[!duplicated(am$uniq.var), ]
        rownames(am) <- am$uniq.var
        vcf.func$AM_pathogenicity <- am[vcf.uniq.var, "am_pathogenicity"]
        vcf.func$AM_class <- am[vcf.uniq.var, "am_class"]
    }

    # add sample GT
    vcf.func.gt <- cbind(vcf.func, rowSums(gt.mat.simple >0, na.rm=T))
    colnames(vcf.func.gt)[ncol(vcf.func.gt)] <- 'samples.variant'
    # add number of samples with non-NA genotype here
    vcf.func.gt <- cbind(vcf.func.gt, rowSums(!is.na(gt.mat.simple), na.rm=T))
    colnames(vcf.func.gt)[ncol(vcf.func.gt)] <- 'samples.genotyped'
    # add in gt matrix
    vcf.func.gt <- cbind(vcf.func.gt, gt.mat.simple)
    # remove rows with no variants after filtering
    vcf.func.gt <- vcf.func.gt[vcf.func.gt$samples.genotyped >0,]
    
    # write out
    write.table(vcf.func.gt, vcf.output.f, sep='\t', quote=F, row.names = F, col.names = T)
}
