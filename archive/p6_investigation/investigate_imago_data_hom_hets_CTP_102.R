# investigate the project CTP_102 data that was parsed with python and loaded up in R
gen_df <- read.table('~/local_data/p6/investigate_sample_data/CTP_102/p6_parsed_haplotypecaller.tsv', 
                     sep='\t', quote='', header=T, row.names = 1, check.names = F)

# stats about the variants 
p6_GGCC_stats <- read.table('~/local_data/p6/p6_GGCC_stats.tsv', sep='\t', quote='', header=T)
p6_GGCC_stats$pos %in% rownames(gen_df)

# full genotyping across all variants!
apply(gen_df, 1, function(x) sum(is.na(x)))

gen_df_p6 <- gen_df[as.character(p6_GGCC_stats[p6_GGCC_stats$type=='p6', "pos"]), ]
# we detect 1,3,4 in GGCC
gen_df_GGCC <- gen_df[as.character(p6_GGCC_stats[p6_GGCC_stats$type=='GGCC', "pos"])[c(1,3,4)], ]
dim(gen_df_p6)
dim(gen_df_GGCC)

# those that are homozygous for p6, what are they for GGCC
p6_hom_test <- c(2,2,2,0,2,2)
p6_hom_test_ind <- which(apply(gen_df_p6,2,function(x) all(x==p6_hom_test)))
gen_df_GGCC[,p6_hom_test_ind]
# most are 222 for GGCC, with some having a zero in the center position in ONE sample, 
# but other samples from the same patient eliminate that, which is strange. 

# remove samples with any NA
gen_df_p6 <- gen_df_p6[, !is.na(colSums(gen_df_p6))]
gen_df_str <- data.frame(
    sample=colnames(gen_df_p6), 
    patient = sapply(colnames(gen_df_p6), function(x) substr(x, 1,7)),
    genotype = apply(gen_df_p6, 2, function(x) paste(as.character(x), collapse=''))
)

# user has alerted me to a few bad patients in this dataset
remove_patients <- c('009-105',
                     '034-101')
gen_df_str <- gen_df_str[!(gen_df_str$patient %in% remove_patients),]

# write a program to phase these results into the most likely constituent haplotypes
# given a list of haplotypes encoded as 0-1 integers and their frequency,
# return the probability of observing the given genotype.
# do this as strings with a layer going to ints for simplicity
# add the new haplotype I see in the project data
ref_df <- data.frame(hap_string = c("000000","111011","011100","111000","011000","111100","000100","100000","001100","000011","011011","010000","010100","001000","110000","100011","001011","111110","110011","101000","101011","100100"),
                    freq = c(0.354,0.2294,0.2222,0.0966,0.0353,0.021,0.0094,0.0082,0.0052,0.0052,0.0048,0.0026,0.0014,0.0012,0.0012,0.0008,0.0004,0.0002,0.0002,0.0002,0.0002,0.0002))

# given two haplotype strings, calculate the integer value from combining them
str_to_genotype <- function(str1, str2){
    ints1 <- as.integer(strsplit(str1, split='')[[1]])
    ints2 <- as.integer(strsplit(str2, split='')[[1]])
    to_return <- paste(as.character(ints1 + ints2), collapse='')
    return(to_return)
}

# calculate this for all pairs
comb_inds <- combn(nrow(ref_df),2)
# need to add doubles
comb_inds <- cbind(comb_inds, rbind(1:nrow(ref_df), 1:nrow(ref_df)))

comb_df <- data.frame(
    hap1 = ref_df[comb_inds[1,],"hap_string"],
    hap2 = ref_df[comb_inds[2,],"hap_string"],
    genotype = apply(comb_inds, 2, function(x){
    str_to_genotype(ref_df[x[1], 'hap_string'],
                    ref_df[x[2], 'hap_string'])
        }
    ),
    exp_freq = apply(comb_inds, 2, function(x){
        ref_df[x[1], 'freq'] * ref_df[x[2], 'freq']
        }
    )
)

# for each sample, what is the most likely haplotype combination
get_haplotype_combinations <- function(g){
    to_return <- comb_df[which(comb_df$genotype==g), ]
    if(nrow(to_return) < 1){
        return(data.frame(hap1=NA,
                          hap2=NA, 
                          genotype=g, 
                          exp_freq=NA, 
                          prob=NA))
    } else {
        # calc probability
        to_return$prob <- round(to_return$exp_freq / sum(to_return$exp_freq),4)
        return(to_return)
    }
}

phase_res <- do.call(rbind, lapply(gen_df_str$genotype, function(x) get_haplotype_combinations(x)[1,]))
phase_res$sample <- gen_df_str$sample
phase_res$patient <- gen_df_str$patient
print(paste('percent of phased samples:',  round(100- (sum(is.na(phase_res$hap1)) / nrow(phase_res)) *100,2)))
hist(phase_res$prob, xlim = c(0,1))
phase_res[order(phase_res$prob)[1:10], ]
sort(table(phase_res[is.na(phase_res$hap1), "genotype"]), decreasing = T)

# take a threshold of 99% probability
phase_prob_thresh <- 0.99
confident_phase_n <- sum(!is.na(phase_res$hap1) & 
                             phase_res$prob >= phase_prob_thresh)
print(paste('percent of confidently phased samples:',  round(confident_phase_n / nrow(phase_res) *100 ,2)))
phase_res_filt <- phase_res[!is.na(phase_res$hap1) & 
                            phase_res$prob >= phase_prob_thresh, ]

# are all samples from one patient phased the same
# this calc done before limiting to one sample per pt
# NO, there could be LOH or other mutations playing into this
for (pt in unique(phase_res_filt$patient)){
    a <- phase_res_filt[phase_res_filt$patient==pt, "genotype"]
    if(!(all(a==a[1]))){
        print( phase_res_filt[phase_res_filt$patient==pt,colnames(phase_res_filt) != 'exp_freq'])
    }
}

# 1 patient do not have consistent genotyping across all samples 
# can just remove some of these samples after manual review
remove_samples <- c('006-102_Dx_PB', '006-102_ITPEOW_BM',
                    '009-104_Dx_PB',
                    '009-105_Dx_PG', '009-105_ET_PG', '009-105_Dx_BS',
                    '011-104_ATP2D168_PG', '011-104_ATP3D168_PG', '011-104_ATP4D168_PG', '011-104_ATPD168_PG', '011-104_Dx_PG', '011-104_ITPD168_PG', '011-104_ITPD84_PG',
                    '011-108_Dx_PG', '011-108_ITPD168_PG', '011-108_ITPD84_PG',
                    '015-105_Dx_PG', '015-105_ITPD84_PG',
                    '031-102_Dx_PG', '031-102_ITPD168_PG', '031-102_ITPD84_PG',
                    '040-101_ATP2D168_PG', '040-101_Dx_PG', '040-101_ITPD168_PG', '040-101_ITPD84_PG',
                    '041-101_Dx_PG',
                    '041-102_ATPEOT_PG', '041-102_Dx_PG', '041-102_ITPD168_PG', '041-102_ITPD84_PG',
                    '060-101_Dx_PG',
                    '060-107_ATP1D168_PG', '060-107_EOS_PG', '060-107_ITPD84_PG',
                    '060-119_ITPD84_PG'
                    )
phase_res_filt <- phase_res_filt[!(phase_res_filt$sample %in% remove_samples),]

# since all patients are consistent, we can just take first representative 
# of a given patient as their genotype / haplotype
phase_res_pt <- phase_res_filt[!duplicated(phase_res_filt$patient), ]
phase_res_pt <- phase_res_pt[order(phase_res_pt$patient), ]

# which of thse have a p6
p6_pattern <- '111011'
p6_hom_pattern <- '222022'
phase_res_pt$p6_positive <- FALSE
phase_res_pt$p6_het <- FALSE
phase_res_pt$p6_hom <- FALSE
# any p6
phase_res_pt[phase_res_pt$hap1==p6_pattern | phase_res_pt$hap2==p6_pattern, 'p6_positive'] <- TRUE
# p6 hom 
phase_res_pt[phase_res_pt$genotype==p6_hom_pattern, 'p6_hom'] <- TRUE
# p6 het 
phase_res_pt[phase_res_pt$p6_positive & phase_res_pt$genotype!=p6_hom_pattern, 'p6_het'] <- TRUE

table(phase_res_pt$p6_positive)
table(phase_res_pt$p6_hom)
table(phase_res_pt$p6_het)

# read the metadata and do the association
metadata <- read.table('~/local_data/p6/investigate_sample_data/CTP_102/CTP_102_metadata.tsv', sep='\t', quote='', header=T, comment.char = '')
rownames(metadata) <- metadata$pt
rownames(phase_res_pt) <- phase_res_pt$patient
table(metadata$JAK2_V617F, useNA = 'ifany')
table(metadata$chr9p_loh, useNA = 'ifany')
table(metadata$any_loh, useNA = 'ifany')
table(metadata$any_chr_event, useNA = 'ifany')

metadata$p6_positive <- phase_res_pt[metadata$pt, "p6_positive"]
metadata$p6_hom <- phase_res_pt[metadata$pt, "p6_hom"]
metadata$p6_het <- phase_res_pt[metadata$pt, "p6_het"]

table(metadata$p6_positive, metadata$JAK2_V617F)
fisher.test(table(metadata$p6_positive, metadata$JAK2_V617F))$p.value
table(metadata$p6_positive, metadata$chr9p_loh)
fisher.test(table(metadata$p6_positive, metadata$chr9p_loh))$p.value
table(metadata$p6_positive, metadata$any_chr_event)
fisher.test(table(metadata$p6_positive, metadata$any_chr_event))$p.value

fisher.test(table(metadata$p6_hom, metadata$JAK2_V617F))$p.value
table(metadata$p6_hom, metadata$chr9p_loh)
fisher.test(table(metadata$p6_hom, metadata$chr9p_loh))$p.value
fisher.test(table(metadata$p6_hom, metadata$any_chr_event))$p.value

fisher.test(table(metadata$p6_het, metadata$JAK2_V617F))$p.value
fisher.test(table(metadata$p6_het, metadata$chr9p_loh))$p.value
fisher.test(table(metadata$p6_het, metadata$any_chr_event))$p.value

sum(is.na(metadata$chr9p_loh))

write.table(metadata, '~/local_data/p6/investigate_sample_data/CTP_102/CTP_102_metadata_p6.tsv', sep='\t', quote=F, row.names = F, col.names = T)
write.table(metadata, '~/cloud_storage/shared_data/shared_p6/CTP_102_metadata_p6.tsv', sep='\t', quote=F, row.names = F, col.names = T)
