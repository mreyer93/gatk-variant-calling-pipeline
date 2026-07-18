# I don't know how to get the phased data for p6 other than this
# it's probably not right, but it will work

p6_vcf_f <- gzfile('~/data/1KG_data/individual_chromosomes/1KG_p6.vcf.gz')
p6 <- read.table(p6_vcf_f, sep='\t',quote='')
# subset to sample columns
p6 <- p6[, 10:ncol(p6)]
# get sample names
colnames(p6) <- read.table('~/data/1KG_data/individual_chromosomes/1KG_p6.fam', sep='\t', quote='', header=F)[,2]

# get separate allele matrices
p6_a1 <- apply(p6, 1:2, function(x) as.numeric(strsplit(x, '|')[[1]][1]))
p6_a2 <- apply(p6, 1:2, function(x) as.numeric(strsplit(x, '|')[[1]][3]))

# concat to str
p6_a1_str <- apply(p6_a1, 2, function(x) paste0(x, collapse=''))
p6_a2_str <- apply(p6_a2, 2, function(x) paste0(x, collapse=''))

# search for the right vars
target_var_str <- '111011'

p6_a1_bool <- p6_a1_str == target_var_str
p6_a2_bool <- p6_a2_str == target_var_str
p6_hom_bool <- p6_a1_bool & p6_a2_bool
p6_het_bool <- p6_a1_bool | p6_a2_bool

sum(p6_het_bool)
sum(p6_hom_bool)
names(which(p6_het_bool))
