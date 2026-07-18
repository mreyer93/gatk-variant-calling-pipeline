# evaluation of the p6 heterozygotes in the EGA data
# in addition to verifying the homozygous patients
library(ggplot2)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggExtra)

# get the EGA sample annotations from the metadata file
ega_metadata <- read.table('~/local_data/EGA_data/ega_metadata_join_roh.tsv', sep='\t', header = T, comment.char = '', quote='')
rownames(ega_metadata) <- ega_metadata$sample
ega_loh <- read.table('~/local_data/EGA_data/ega_loh_metadata.tsv', sep='\t', header = T, comment.char = '', quote='')
rownames(ega_loh) <- ega_loh$patient
# metadata about if the PATIENT has JAK2_V617F
ega_metadata_patient <- ega_metadata[!duplicated(ega_metadata$patient),]
ega_metadata_patient$JAK2_V617F <- sapply(ega_metadata_patient$patient, function(x) any(ega_metadata[ega_metadata$patient==x, "has_jak2_V617F"]))
rownames(ega_metadata_patient) <- ega_metadata_patient$patient

# genotype matrix
genotype <- read.table('~/local_data/p6/p6_pattern/genotype_matrix.tsv', sep='\t', quote='', row.names = 1, header=T, check.names = F)

# only interested in the normal samples right now
genotype_normal <- genotype[, ega_metadata[ega_metadata$tumor_normal=='N', "sample"]]
# or what if we include one sample from each patient, either the normal or tumor if it's unavaiable
sub_samples <- c(ega_metadata[ega_metadata$tumor_normal=='N', 'sample'],
                 ega_metadata[ega_metadata$tumor.normal.pair == "False", 'sample'])
# one duplicate here
sub_samples <- sub_samples[sub_samples != "PD6647a2"]
genotype_single <- genotype[, sub_samples]

# convert to string
gen_df_str <- data.frame(
    sample=colnames(genotype_single), 
    patient = sapply(colnames(genotype_single), function(x) ega_metadata[x, "patient"]),
    genotype = apply(genotype_single, 2, function(x) paste(as.character(x), collapse=''))
)
# there's a few duplicate patients
dp <- gen_df_str[duplicated(gen_df_str$patient),"patient"]
gen_df_str[gen_df_str$patient %in% dp, ]
gen_df_str <- gen_df_str[!duplicated(gen_df_str$patient),]

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

# since only one sample per patient, don't need to subset here
phase_res_pt <- phase_res_filt

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

# by disease
phase_res_pt$study <- ega_metadata[phase_res_pt$sample, "study"]
table(phase_res_pt$p6_positive, phase_res_pt$study)[2,] / colSums(table(phase_res_pt$p6_positive, phase_res_pt$study))
table(phase_res_pt$p6_hom, phase_res_pt$study)[2,] / colSums(table(phase_res_pt$p6_positive, phase_res_pt$study))
table(phase_res_pt$p6_het, phase_res_pt$study)[2,] / colSums(table(phase_res_pt$p6_positive, phase_res_pt$study))

# add loh information
phase_res_pt$JAK2_V617F <- ega_metadata_patient[phase_res_pt$patient, "JAK2_V617F"]
phase_res_pt$chr9p_loh <- ega_loh[phase_res_pt$patient, "chr9p_loh"]
phase_res_pt$any_loh <- ega_loh[phase_res_pt$patient, "any_loh"]
phase_res_pt$any_chr_event <- ega_loh[phase_res_pt$patient, "any_chr_event"]

# patients in original set that have the T/N pair
pt_with_normal <- intersect(unique(ega_metadata$patient[ega_metadata$tumor.normal.pair =='True']), 
                            phase_res_pt$patient)
rownames(phase_res_pt) <- phase_res_pt$patient

# do FET association
table(phase_res_pt$p6_positive, phase_res_pt$JAK2_V617F, useNA = 'ifany')
fisher.test(table(phase_res_pt$p6_positive, phase_res_pt$JAK2_V617F))$p.value
table(phase_res_pt$p6_positive, phase_res_pt$chr9p_loh, useNA = 'ifany')
fisher.test(table(phase_res_pt$p6_positive, phase_res_pt$chr9p_loh))$p.value
table(phase_res_pt$p6_positive, phase_res_pt$any_loh, useNA = 'ifany')
fisher.test(table(phase_res_pt$p6_positive, phase_res_pt$any_loh))$p.value
table(phase_res_pt$p6_positive, phase_res_pt$any_chr_event, useNA = 'ifany')
fisher.test(table(phase_res_pt$p6_positive, phase_res_pt$any_chr_event))$p.value

p_df <- data.frame(all=c(
    fisher.test(table(phase_res_pt$p6_positive, phase_res_pt$JAK2_V617F))$p.value,
    fisher.test(table(phase_res_pt$p6_positive, phase_res_pt$chr9p_loh))$p.value,
    fisher.test(table(phase_res_pt$p6_positive, phase_res_pt$any_loh))$p.value,
    fisher.test(table(phase_res_pt$p6_positive, phase_res_pt$any_chr_event))$p.value,
    
    fisher.test(table(phase_res_pt$p6_hom, phase_res_pt$JAK2_V617F))$p.value,
    fisher.test(table(phase_res_pt$p6_hom, phase_res_pt$chr9p_loh))$p.value,
    fisher.test(table(phase_res_pt$p6_hom, phase_res_pt$any_loh))$p.value,
    fisher.test(table(phase_res_pt$p6_hom, phase_res_pt$any_chr_event))$p.value,
    
    fisher.test(table(phase_res_pt$p6_het, phase_res_pt$JAK2_V617F))$p.value,
    fisher.test(table(phase_res_pt$p6_het, phase_res_pt$chr9p_loh))$p.value,
    fisher.test(table(phase_res_pt$p6_het, phase_res_pt$any_loh))$p.value,
    fisher.test(table(phase_res_pt$p6_het, phase_res_pt$any_chr_event))$p.value
), 
    ET=c(
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_positive'], phase_res_pt[phase_res_pt$study=='ET', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_positive'], phase_res_pt[phase_res_pt$study=='ET', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_positive'], phase_res_pt[phase_res_pt$study=='ET', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_positive'], phase_res_pt[phase_res_pt$study=='ET', 'any_chr_event']))$p.value,
    
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_hom'], phase_res_pt[phase_res_pt$study=='ET', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_hom'], phase_res_pt[phase_res_pt$study=='ET', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_hom'], phase_res_pt[phase_res_pt$study=='ET', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_hom'], phase_res_pt[phase_res_pt$study=='ET', 'any_chr_event']))$p.value,
    
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_het'], phase_res_pt[phase_res_pt$study=='ET', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_het'], phase_res_pt[phase_res_pt$study=='ET', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_het'], phase_res_pt[phase_res_pt$study=='ET', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='ET', 'p6_het'], phase_res_pt[phase_res_pt$study=='ET', 'any_chr_event']))$p.value
),
    PMF=c(
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_positive'], phase_res_pt[phase_res_pt$study=='PMF', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_positive'], phase_res_pt[phase_res_pt$study=='PMF', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_positive'], phase_res_pt[phase_res_pt$study=='PMF', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_positive'], phase_res_pt[phase_res_pt$study=='PMF', 'any_chr_event']))$p.value,
    
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_hom'], phase_res_pt[phase_res_pt$study=='PMF', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_hom'], phase_res_pt[phase_res_pt$study=='PMF', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_hom'], phase_res_pt[phase_res_pt$study=='PMF', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_hom'], phase_res_pt[phase_res_pt$study=='PMF', 'any_chr_event']))$p.value,
    
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_het'], phase_res_pt[phase_res_pt$study=='PMF', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_het'], phase_res_pt[phase_res_pt$study=='PMF', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_het'], phase_res_pt[phase_res_pt$study=='PMF', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PMF', 'p6_het'], phase_res_pt[phase_res_pt$study=='PMF', 'any_chr_event']))$p.value
),
    PV=c(
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_positive'], phase_res_pt[phase_res_pt$study=='PV', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_positive'], phase_res_pt[phase_res_pt$study=='PV', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_positive'], phase_res_pt[phase_res_pt$study=='PV', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_positive'], phase_res_pt[phase_res_pt$study=='PV', 'any_chr_event']))$p.value,
    
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_hom'], phase_res_pt[phase_res_pt$study=='PV', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_hom'], phase_res_pt[phase_res_pt$study=='PV', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_hom'], phase_res_pt[phase_res_pt$study=='PV', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_hom'], phase_res_pt[phase_res_pt$study=='PV', 'any_chr_event']))$p.value,
    
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_het'], phase_res_pt[phase_res_pt$study=='PV', 'JAK2_V617F']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_het'], phase_res_pt[phase_res_pt$study=='PV', 'chr9p_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_het'], phase_res_pt[phase_res_pt$study=='PV', 'any_loh']))$p.value,
    fisher.test(table(phase_res_pt[phase_res_pt$study=='PV', 'p6_het'], phase_res_pt[phase_res_pt$study=='PV', 'any_chr_event']))$p.value
)
)

rownames(p_df) <- c('P6 positive, JAK2','P6 positive, chr9p_LOH','P6 positive, any_LOH','P6 positive, any_chr_event','P6 hom, JAK2','P6 hom, chr9p_LOH','P6 hom, any_LOH','P6 hom, any_chr_event','P6 het, JAK2','P6 het, chr9p_LOH','P6 het, any_LOH','P6 het, any_chr_event')
p_df


# cases of homozygotes
hom_match = c(2,2,2,0,2,2)
names(hom_match) <- rownames(genotype_normal)
hom_index <- which(apply(genotype_normal, 2, function(x) identical(x, hom_match)))
print(length(hom_index))
# cases of homozygotes in genotype_single
hom_index_single <- which(apply(genotype_single, 2, function(x) identical(x, hom_match)))
print(length(hom_index_single))

# homozygotes of the opposite allele
hom_opp_match = c(0,0,0,2,0,0)
names(hom_opp_match) <- rownames(genotype_normal)
hom_opp_index <- which(apply(genotype_normal, 2, function(x) identical(x, hom_opp_match)))
print(length(hom_opp_index))

# heteozygotes for everything
het_match = c(1,1,1,1,1,1)
names(het_match) <- rownames(genotype_normal)
het_index <- which(apply(genotype_normal, 2, function(x) identical(x, het_match)))
print(length(het_index))

# most popular haplotype in Europeans is all reference
pop_hom_match <- c(0,0,0,0,0,0)
names(pop_hom_match) <- rownames(genotype_normal)
print(length(which(apply(genotype_normal, 2, function(x) identical(x, pop_hom_match)))))
# as a p6 het with this you have mostly hets. Except 4th allele which is always ref
pop_het_match <- c(1,1,1,0,1,1)
names(pop_het_match) <- rownames(genotype_normal)
print(length(which(apply(genotype_normal, 2, function(x) identical(x, pop_het_match)))))

# count presence of any p6 haplotype
# it can be one of 4 things I can detect in the data
# (based on my google sheet tracking this)
p6_possibilities <- list(c(1,1,1,0,1,1),
                         c(2,2,2,0,2,2),
                         c(1,2,2,1,1,1),
                         c(2,1,1,0,1,1))

# need a T/F if each patient is positive for ANY of these
genotype_normal_nn <- genotype_normal
rownames(genotype_normal_nn) <- NULL

p6_haplotype_any <- apply(genotype_normal_nn, 2, function(x) {
    any(
        sapply(p6_possibilities, function(y) {
                identical(x,y)
            }
        )
    )
})

ega_metadata_normal <- ega_metadata[ega_metadata$tumor_normal=='N', ]
ega_metadata_normal$p6_haplotype_any <- p6_haplotype_any[ega_metadata_normal$sample]
ega_metadata_normal$chr9p_loh <- ega_loh[ega_metadata_normal$patient, "chr9p_loh"]
ega_metadata_normal$any_loh <- ega_loh[ega_metadata_normal$patient, "any_loh"]
ega_metadata_normal$any_chr_event <- ega_loh[ega_metadata_normal$patient, "any_chr_event"]

print(nrow(ega_metadata_normal))
table(ega_metadata_normal$p6_haplotype_any, ega_metadata_normal$has_jak2_V617F)
fisher.test(table(ega_metadata_normal$p6_haplotype_any, ega_metadata_normal$has_jak2_V617F))

table(ega_metadata_normal$p6_haplotype_any, ega_metadata_normal$chr9p_loh)
fisher.test(table(ega_metadata_normal$p6_haplotype_any, ega_metadata_normal$chr9p_loh))

table(ega_metadata_normal$p6_haplotype_any, ega_metadata_normal$any_loh)
fisher.test(table(ega_metadata_normal$p6_haplotype_any, ega_metadata_normal$any_loh))

table(ega_metadata_normal$p6_haplotype_any, ega_metadata_normal$any_chr_event)
fisher.test(table(ega_metadata_normal$p6_haplotype_any, ega_metadata_normal$any_chr_event))

# can we get frequency of columns
gl <- apply(genotype_normal, 2, function(x) paste(x, collapse='_'))
gt <- sort(table(gl), decreasing = T)
gdf <- data.frame(gt)
gdf$pct = gdf$Freq / sum(gdf$Freq)

############################################################
# SPLIT BY STUDY
############################################################
# MF
############################################################
ega_metadata_normal_PMF <- ega_metadata_normal[ega_metadata$study=='PMF', ]
ega_metadata_normal_PMF <- ega_metadata_normal_PMF[!is.na(ega_metadata_normal_PMF$study),]
print(nrow(ega_metadata_normal_PMF))
table(ega_metadata_normal_PMF$p6_haplotype_any, ega_metadata_normal_PMF$has_jak2_V617F)
fisher.test(table(ega_metadata_normal_PMF$p6_haplotype_any, ega_metadata_normal_PMF$has_jak2_V617F))
table(ega_metadata_normal_PMF$p6_haplotype_any, ega_metadata_normal_PMF$chr9p_loh)
fisher.test(table(ega_metadata_normal_PMF$p6_haplotype_any, ega_metadata_normal_PMF$chr9p_loh))
table(ega_metadata_normal_PMF$p6_haplotype_any, ega_metadata_normal_PMF$any_loh)
fisher.test(table(ega_metadata_normal_PMF$p6_haplotype_any, ega_metadata_normal_PMF$any_loh))
table(ega_metadata_normal_PMF$p6_haplotype_any, ega_metadata_normal_PMF$any_chr_event)
fisher.test(table(ega_metadata_normal_PMF$p6_haplotype_any, ega_metadata_normal_PMF$any_chr_event))
############################################################
# ET
############################################################
ega_metadata_normal_ET <- ega_metadata_normal[ega_metadata$study=='ET', ]
ega_metadata_normal_ET <- ega_metadata_normal_ET[!is.na(ega_metadata_normal_ET$study),]
print(nrow(ega_metadata_normal_ET))
table(ega_metadata_normal_ET$p6_haplotype_any, ega_metadata_normal_ET$has_jak2_V617F)
fisher.test(table(ega_metadata_normal_ET$p6_haplotype_any, ega_metadata_normal_ET$has_jak2_V617F))
table(ega_metadata_normal_ET$p6_haplotype_any, ega_metadata_normal_ET$chr9p_loh)
fisher.test(table(ega_metadata_normal_ET$p6_haplotype_any, ega_metadata_normal_ET$chr9p_loh))
table(ega_metadata_normal_ET$p6_haplotype_any, ega_metadata_normal_ET$any_loh)
fisher.test(table(ega_metadata_normal_ET$p6_haplotype_any, ega_metadata_normal_ET$any_loh))
table(ega_metadata_normal_ET$p6_haplotype_any, ega_metadata_normal_ET$any_chr_event)
fisher.test(table(ega_metadata_normal_ET$p6_haplotype_any, ega_metadata_normal_ET$any_chr_event))
############################################################
# PV
############################################################
ega_metadata_normal_PV <- ega_metadata_normal[ega_metadata$study=='PV', ]
ega_metadata_normal_PV <- ega_metadata_normal_PV[!is.na(ega_metadata_normal_PV$study),]
print(nrow(ega_metadata_normal_PV))
table(ega_metadata_normal_PV$p6_haplotype_any, ega_metadata_normal_PV$has_jak2_V617F)
fisher.test(table(ega_metadata_normal_PV$p6_haplotype_any, ega_metadata_normal_PV$has_jak2_V617F))
table(ega_metadata_normal_PV$p6_haplotype_any, ega_metadata_normal_PV$chr9p_loh)
fisher.test(table(ega_metadata_normal_PV$p6_haplotype_any, ega_metadata_normal_PV$chr9p_loh))
table(ega_metadata_normal_PV$p6_haplotype_any, ega_metadata_normal_PV$any_loh)
fisher.test(table(ega_metadata_normal_PV$p6_haplotype_any, ega_metadata_normal_PV$any_loh))
table(ega_metadata_normal_PV$p6_haplotype_any, ega_metadata_normal_PV$any_chr_event)
fisher.test(table(ega_metadata_normal_PV$p6_haplotype_any, ega_metadata_normal_PV$any_chr_event))
