# can we read a VCF file and predict the LOH regions from it directly
# without doing a normal comparison
library(vcfR)
library(ggplot2)

# work with the patient 003-203
# a simple case of two samples, where the tumor should have 9P GC
# 003-203_Dx_BS 003-203_Dx_PG
f_normal <- '~/data/CTP_201/02_variants_reference/01_gatk_variant_calling/filtered/003-203_Dx_BS_filtered.vcf'
f_tumor <- '~/data/CTP_201/02_variants_reference/01_gatk_variant_calling/filtered/003-203_Dx_PG_filtered.vcf'

vcf_n <- read.vcfR(f_tumor, verbose = FALSE )
vcf_t <- read.vcfR(f_tumor, verbose = FALSE )

pass_filter <- vcf_n@fix[,'FILTER'] == 'PASS'
valid_gt <- c('0/1', '0|1', '1|0', '1/0')
gt_n <- extract.gt(vcf_n)
keep_mask <- (gt_n %in% valid_gt) & pass_filter
keep_names <- rownames(gt_n)[gt_n %in% valid_gt & pass_filter]
length(keep_names)

af_n <- extract.gt(vcf_n, element = 'AF', as.numeric = T)
dp_n <- extract.gt(vcf_n, element = 'DP', as.numeric = T)

gt_n <- gt_n[keep_names,,drop=F]
af_n <- af_n[keep_names,,drop=F]
dp_n <- dp_n[keep_names,,drop=F]

plot_df <- data.frame(vcf_n@fix[, c('CHROM', 'POS')])
plot_df$POS <- as.numeric(plot_df$POS)
plot_df <- plot_df[keep_mask, ]
plot_df$GT <- gt_n[,1]
plot_df$AF <- af_n[,1]
plot_df$DP <- dp_n[,1]

min_dp <- 10
plot_df <- plot_df[plot_df$DP >= min_dp, ]
plot_df$i <- 1:nrow(plot_df)

ggplot(plot_df, aes(x=i, y=AF)) + 
    geom_point()

ggplot(plot_df[plot_df$CHROM=='chr9',], aes(x=POS, y=AF)) + 
    geom_point()
