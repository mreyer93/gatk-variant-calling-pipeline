# investigate SF3B1 in the project data
# for the three trials

# load processed germline data
f.101 <- '~/local_data/CTP_101/07_joint_vcf/02_variant_annotations/annotations_combined.vcf.gz'
f.102 <- '~/local_data/CTP_102/07_joint_vcf/02_variant_annotations/annotations_combined.vcf.gz'
f.201 <- '~/local_data/CTP_201/07_joint_vcf/02_variant_annotations/annotations_combined.vcf.gz'

dat.101 <- read.table(gzfile(f.101), sep='\t', quote='', header=T, check.names = F)
dat.102 <- read.table(gzfile(f.102), sep='\t', quote='', header=T, check.names = F)
dat.201 <- read.table(gzfile(f.201), sep='\t', quote='', header=T, check.names = F)
dat.list <- list(dat.101, dat.102, dat.201)

# subset to this gene list
gene.list <- c('SF3B1')
dat.list.filter <- lapply(dat.list, function(x) x[x$hugoSymbol %in% gene.list, ])

View(dat.list.filter[[1]])
View(dat.list.filter[[2]])
View(dat.list.filter[[3]])

# not to much here. No common nonsynonmyous variants of note.