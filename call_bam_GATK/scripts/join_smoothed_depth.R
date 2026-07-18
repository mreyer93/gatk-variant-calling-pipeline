#!/usr/bin/env Rscript
# join smoothed depth files

depth_files <- unlist(snakemake@input)
sample_list <- snakemake@params[['sample_list']]
out_f <- snakemake@output[[1]]

r1 <- read.table(depth_files[1], sep='\t', quote='', header=T)
r1 <- r1[, colnames(r1) != 'max_w_depth']

max_w_depth_mat <- do.call(cbind, lapply(depth_files, function(f){
    read.table(f, sep='\t', quote='', header=T)[, 'max_w_depth']
}))

print(dim(max_w_depth_mat))
colnames(max_w_depth_mat) <- sample_list

out_df <- cbind(r1, max_w_depth_mat)
write.table(out_df, out_f, sep='\t', row.names=F, col.names=T, quote=F)