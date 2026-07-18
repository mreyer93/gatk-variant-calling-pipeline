# investigate the project CTP_201 data that was parsed with python and loaded up in R
gen_df <- read.table('~/local_data/p6/investigate_sample_data/CTP_201/p6_parsed_haplotypecaller.tsv', 
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

# only keep the data relating to one NORMAL sample per patient, or the best sample we have
# keep.samples <- c("001-201_Dx_BS","001-203_Dx_PG","001-204_Dx_HS","001-205_Dx_BS","001-206_Dx_PG","001-207_Dx_HS","003-201_Dx_BS","003-202_Dx_HS","003-203_Dx_BS","003-204_Dx_BS","003-205_Dx_HS","003-206_Dx_BS","006-201_Dx_HS","006-202_Dx_BS","008-201_Dx_HS","009-201_Dx_BS","009-202_Dx_BS","012-201_Dx_BS","012-202_Dx_BS","014-201_Dx_HS","014-202_Dx_HS","014-203_Dx_PG","014-204_Dx_HS","014-205_Dx_HS","014-206_Dx_HS","020-201_Dx_HS","020-202_Dx_HS","020-203_Dx_HS","020-204_Dx_PG","020-205_Dx_HS","021-201_Dx_HS","024-201_Dx_BS","024-202_Dx_HS","024-203_Dx_HS","024-204_Dx_HS","024-205_Dx_PG","030-201_Dx_BS","030-202_Dx_BS","031-201_Dx_HS","031-202_Dx_HS","032-201_Dx_BS","032-202_Dx_HS","032-203_Dx_BS","032-204_Dx_BS","032-205_Dx_BS","032-206_Dx_BS","034-201_Dx_BS","034-202_Dx_HS","034-203_Dx_PG","034-204_Dx_BS","034-205_Dx_PG","034-206_Dx_HS","040-201_Dx_HS","040-203_Dx_HS","040-204_Dx_HS","050-201_Dx_HS","050-202_Dx_BS","050-204_Dx_HS","050-205_Dx_HS","051-201_Dx_BS","051-202_Dx_HS","051-203_Dx_HS","055-201_Dx_BS","055-202_Dx_BS","055-203_Dx_BS","055-204_Dx_BS","055-205_Dx_BS","055-206_Dx_BS","055-207_Dx_HS","056-201_Dx_BS","056-202_Dx_HS","057-201_Dx_HS","060-201_Dx_HS","060-202_Dx_HS","060-203_Dx_HS","060-204_Dx_HS","060-205_Dx_HS","060-206_Dx_HS","060-207_Dx_HS","060-208_Dx_HS","060-209_Dx_BS","060-210_Dx_HS","060-211_Dx_PG","060-212_Dx_HS","060-213_Dx_HS","060-214_Dx_HS","060-215_Dx_HS","060-216_Dx_HS","060-217_Dx_PG","060-218_Dx_HS")
# all(keep.samples %in% colnames(gen_df))
# gen_df <- gen_df[, keep.samples]
# dim(gen_df)
# gen_df <- read.table('~/local_data/p6/p6_pattern/genotype_matrix.tsv', 
#                      sep='\t', quote='', header=T, row.names = 1, check.names = F)

# separate into p6 and GGCC when necessary

# remove samples with any NA
gen_df_p6 <- gen_df_p6[, !is.na(colSums(gen_df_p6))]
gen_df_str <- data.frame(
    sample=colnames(gen_df_p6), 
    patient = sapply(colnames(gen_df_p6), function(x) substr(x, 1,7)),
    genotype = apply(gen_df_p6, 2, function(x) paste(as.character(x), collapse=''))
)

# user has alerted me to a few bad patients in this dataset
remove_patients <- c('001-202',
                     '050-203',
                     '056-203')
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
        print( phase_res_filt[phase_res_filt$patient==pt,])
    }
}

# 4 patients do not have consistent genotyping across all samples 
# can just remove some of these samples after manual review
remove_samples <- c('055-201_Dx_HS')
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

table(phase_res_pt$p6_positive, useNA = 'ifany')
table(phase_res_pt$p6_hom, useNA = 'ifany')
table(phase_res_pt$p6_het, useNA = 'ifany')

# read the metadata and do the association
metadata <- read.table('~/local_data/p6/investigate_sample_data/CTP_201/CTP_201_metadata.tsv', sep='\t', quote='', header=T, comment.char = '')
rownames(metadata) <- metadata$pt
rownames(phase_res_pt) <- phase_res_pt$patient
table(metadata$JAK2_V617F, useNA = 'ifany')
table(metadata$chr9p_loh, useNA = 'ifany')
table(metadata$any_loh, useNA = 'ifany')
table(metadata$any_chr_event, useNA = 'ifany')

metadata$p6_positive <- phase_res_pt[metadata$pt, "p6_positive"]
metadata$p6_hom <- phase_res_pt[metadata$pt, "p6_hom"]
metadata$p6_het <- phase_res_pt[metadata$pt, "p6_het"]

# Test for associations with P6 positive
table(metadata$p6_positive, metadata$JAK2_V617F)
fisher.test(table(metadata$p6_positive, metadata$JAK2_V617F))$p.value
table(metadata$p6_positive, metadata$chr9p_loh, useNA = 'ifany')
fisher.test(table(metadata$p6_positive, metadata$chr9p_loh))$p.value
table(metadata$p6_positive, metadata$any_chr_event)
fisher.test(table(metadata$p6_positive, metadata$any_chr_event))$p.value

# Test for associations with P6 homozygous
fisher.test(table(metadata$p6_hom, metadata$JAK2_V617F))$p.value
table(metadata$p6_hom, metadata$chr9p_loh, useNA = 'ifany')
fisher.test(table(metadata$p6_hom, metadata$chr9p_loh))$p.value
fisher.test(table(metadata$p6_hom, metadata$any_chr_event))$p.value

# Test for associations with P6 heterozygous
fisher.test(table(metadata$p6_het, metadata$JAK2_V617F))$p.value
fisher.test(table(metadata$p6_het, metadata$chr9p_loh))$p.value
fisher.test(table(metadata$p6_het, metadata$any_chr_event))$p.value

table(metadata$JAK2_V617F, metadata$chr9p_loh, useNA = 'ifany')
fisher.test(table(metadata$JAK2_V617F, metadata$chr9p_loh))

# write out the results of the phasing
write.table(metadata, '~/local_data/p6/investigate_sample_data/CTP_201/CTP_201_metadata_p6.tsv', sep='\t', quote=F, row.names = F, col.names = T)
write.table(metadata, '~/cloud_storage/shared_data/shared_p6/CTP_201_metadata_p6.tsv', sep='\t', quote=F, row.names = F, col.names = T)
