# investigate JAK2 and GFI1 variants in the EGA data
# the files are so large that the automatic combination of anntation files 
# did not complete, so I have grepped for the genes and will combine the 
# files manually here. 
depth.min <- 5

vcf.header.f <- '~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/funcotator_header.vcf'
vcf.func.GFI1.f <- '~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/funcotator_GFI1.vcf'
vcf.func.JAK2.f <- '~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/funcotator_JAK2.vcf'
vcf.cadd.header.f <- '~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/CADD_simple_table_header.vcf'
vcf.cadd.GFI1.f <- '~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/CADD_simple_table_GFI1.vcf'
vcf.cadd.JAK2.f <- '~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/CADD_simple_table_JAK2.vcf'
vcf.cadd.full.GFI1.f <- '~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/CADD_full_table_GFI1.vcf'
vcf.cadd.full.JAK2.f <- '~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/CADD_full_table_JAK2.vcf'

vcf.header <- read.table(vcf.header.f, sep='\t', quote='', comment.char = '', header = F, fill=T)
vcf.header[1,1] <- 'CHROM'

vcf.func.GFI1 <- read.table(vcf.func.GFI1.f, sep='\t', quote='', comment.char = '', header = F)
vcf.func.JAK2 <- read.table(vcf.func.JAK2.f, sep='\t', quote='', comment.char = '', header = F)
vcf <- rbind(vcf.func.GFI1, vcf.func.JAK2)
colnames(vcf) <- vcf.header[1,]

vcf.cadd.GFI1 <- read.table(vcf.cadd.GFI1.f, sep='\t', quote='', comment.char = '', header = F)
vcf.cadd.JAK2 <- read.table(vcf.cadd.JAK2.f, sep='\t', quote='', comment.char = '', header = F)
vcf.cadd.full.GFI1 <- read.table(vcf.cadd.full.GFI1.f, sep='\t', quote='', comment.char = '', header = F)
vcf.cadd.full.JAK2 <- read.table(vcf.cadd.full.JAK2.f, sep='\t', quote='', comment.char = '', header = F)
cadd <- rbind(vcf.cadd.GFI1, vcf.cadd.JAK2)
cadd.full <- rbind(vcf.cadd.full.GFI1, vcf.cadd.full.JAK2)
colnames(cadd) <- c("CHROM", "Pos", "Ref", "Alt", "Type", "Length", "ConsDetail", "GeneName", "PHRED")
colnames(cadd.full) <- read.table('~/data/EGA_data_whole_exome/07_joint_vcf/02_variant_annotations/CADD_full_header.txt')[1,]

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

# add CADD information
colnames(cadd)[1] <- 'CHROM'
# unique names for cadd variants and other variants
cadd$uniq.var <- paste(gsub('chr', '', cadd$CHROM), cadd$Pos, cadd$Ref, cadd$Alt, sep='_')
# remove duplicates here because all I care about is the phred score
cadd <- cadd[!duplicated(cadd$uniq.var), ]
rownames(cadd) <- cadd$uniq.var
vcf.uniq.var <- paste(gsub('chr', '', vcf.func$CHROM), vcf.func$POS, vcf.func$REF, vcf.func$ALT, sep='_')
vcf.func$CADD_phred <- cadd[vcf.uniq.var, "PHRED"]

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
vcf.func.gt <- vcf.func.gt[vcf.func.gt$samples.variant >0,]
View(vcf.func.gt)
