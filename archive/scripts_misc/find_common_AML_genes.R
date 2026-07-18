# something to look for interesting variants across a large set of VCFs
# and report the stats I want on them
# takes in metadata...
# reports in long form, and then perhaps a wide form with lots of columns?
# SAMPLE PATIENT TN GENE VARIANT CADD AF DEPTH ETC

metadata.f <- '~/local_data/EGA_data/bam_metadata_uniq_samp.tsv'
annot_combined_dir <- '~/local_data/EGA_data/02_variants_reference/annotations_combined/'
outf <- '~/local_data/EGA_data/02_variants_reference/common_gene_findings.tsv'
min.cadd.report <- 20

# input directories
indir.unfilt <- file.path(annot_combined_dir, 'unfiltered')
# read metadata
metadata <- read.table(metadata.f, sep='\t', header=T, comment.char ='', fill=T, quote='')
rownames(metadata) <- metadata$sample

report.df.big <- data.frame()
for (sample in metadata$sample){
    print(sample)
    vcf.unfilt.f <- file.path(indir.unfilt, paste0(sample, '_unfiltered_annotated.vcf'))
    if(!all(file.exists(vcf.unfilt.f))){
        print(vcf.unfilt.f)
        stop('Some vcf files do not exist')
    }

    vcf.df <- lapply(vcf.unfilt.f, function(x) {
        a <- read.table(x, header=T, sep='\t', quote='', comment.char = '', check.names = F)
        a <- a[order(a$CHROM, a$POS, a$REF, a$ALT),]
        a$key <- paste(a$CHROM, a$POS, a$REF, a$ALT, sep='_')
        return(a)
    })[[1]]

    vcf.df$proteinChange[is.na(vcf.df$proteinChange)] <- ''
    vcf.df$gene_protein <- paste(vcf.df$hugoSymbol, vcf.df$proteinChange, sep='_')
    vcf.df$gene_protein[vcf.df$proteinChange==''] <- paste(vcf.df$hugoSymbol[vcf.df$proteinChange==''],
                                                             vcf.df$variantClassification[vcf.df$proteinChange==''], sep='_')

    # sort by cadd and eliminate those less than threshold
    vcf.df <- vcf.df[vcf.df$CADD_phred >= min.cadd.report, ]
    vcf.df <- vcf.df[order(vcf.df$CADD_phred, decreasing = T), ]

    search.specific <- c('JAK2_p.V617F')
    search.gene <- c('JAK2', 'CALR', 'MPL')

    specific.df <- do.call(rbind, lapply(search.specific, function(x) vcf.df[grep(x, vcf.df$gene_protein), ]))
    gene.df <- do.call(rbind, lapply(search.gene, function(x) vcf.df[grep(x, vcf.df$hugoSymbol), ]))

    report.df <- rbind(specific.df, gene.df)
    report.df <- report.df[!duplicated(report.df$key),]

    if(nrow(report.df)==0){
        print('no variants!')
    } else{

        # make a version with only the right info
        colnames(report.df)[colnames(report.df)==sample] <- 'AF'
        colnames(report.df)[colnames(report.df)=='key'] <- 'genome_change'
        report.df$sample <- sample
        report.df$patient <- metadata[sample, "patient"]
        report.df$tumor_normal <- metadata[sample, "tumor_normal"]

        keep.cols <- c('sample', 'patient', 'tumor_normal', 'genome_change', 'hugoSymbol', 'gene_protein', 'CADD_phred', 'AF')
        report.df.sub <- report.df[, keep.cols]
        report.df.big <- rbind(report.df.big, report.df.sub)
    }

}

report.df.big <- do.call(rbind,report.list)

write.table(report.df.big, outf, sep='\t', quote=F, row.names=F, col.names=T)
