#!/usr/bin/env Rscript
# new LOH calling script on the project data
# working with the somatic variant calls, and depth from a single sample
#
# For each target region: the highest mean depth over 200 bp windows stepped by 100 bp
# (the mean over the whole region when it is shorter than one window). The windows are
# laid out by the depth_windows rule and their depth summed by `samtools bedcov`, so this
# script only takes a maximum per region.
#
# It used to read the per-base `samtools depth -a` output into R and scan a whole
# chromosome's rows for every region, which on an exome (~134 million positions, ~200,000
# regions) never finished. It also set the region names as row names, which fails as soon
# as a name repeats (exon BEDs often name each exon after its gene), and gave 0 to any
# region exactly one window long.

depth_f <- snakemake@input[[1]]
regions_f <- snakemake@input[[2]]
out_f <- snakemake@output[[1]]

# windowed depth: chr, start, end, region (its line number in the target file), and the
# sum of per-base depth over the window
win_df <- read.table(gzfile(depth_f), header=F, sep='\t',
                     colClasses=c('character', 'integer', 'integer', 'integer', 'numeric'),
                     col.names=c('chr', 'start', 'end', 'region', 'depth_sum'))
win_df$mean_depth <- win_df$depth_sum / (win_df$end - win_df$start)
max_w_depth <- tapply(win_df$mean_depth, win_df$region, max)

# read regions file, skipping the same lines depth_windows skips so line numbers match.
# A 3-column BED is accepted too: strand '+' and a chr:start-end name are filled in.
region_lines <- readLines(regions_f)
region_lines <- region_lines[nzchar(trimws(region_lines)) &
                             !grepl('^(#|track|browser)', region_lines)]
regions_df <- read.table(text=region_lines, sep='\t', quote='', header=F,
                         comment.char='', stringsAsFactors=F)
if (ncol(regions_df) < 5) {
    regions_df$V4 <- '+'
    regions_df$V5 <- paste0(regions_df$V1, ':', regions_df$V2, '-', regions_df$V3)
}
regions_df <- regions_df[, 1:5]
colnames(regions_df) <- c('chr', 'start', 'end', 'strand', 'name')
regions_df$max_w_depth <- round(as.numeric(max_w_depth[as.character(seq_len(nrow(regions_df)))]))

# numeric chromosomes in order, then X, Y, M, then anything else in file order
chr_rank <- suppressWarnings(as.numeric(gsub('chr', '', regions_df$chr)))
chr_rank[regions_df$chr %in% c('chrX', 'X')] <- 1001
chr_rank[regions_df$chr %in% c('chrY', 'Y')] <- 1002
chr_rank[regions_df$chr %in% c('chrM', 'chrMT', 'M', 'MT')] <- 1003
other <- is.na(chr_rank)
chr_rank[other] <- 2000 + match(regions_df$chr[other], unique(regions_df$chr[other]))
regions_df <- regions_df[order(chr_rank, regions_df$start), ]
regions_df$chr <- factor(regions_df$chr, levels=unique(regions_df$chr))
regions_df$size <- regions_df$end - regions_df$start
regions_df$i <- 1:nrow(regions_df)
# within-chromosome index of each region. lapply, not sapply: when every chromosome has
# the same number of regions (e.g. exactly one, as with whole-contig targets) sapply
# simplifies to a vector and do.call(c, .) then fails with "second argument must be a list".
regions_df$chr_i <- unlist(lapply(table(regions_df$chr), function(x) seq_len(x)),
                           use.names = FALSE)

smooth_region_df <- regions_df[, c('chr', 'start', 'end', 'strand', 'name', 'size', 'i',
                                   'chr_i', 'max_w_depth')]
# clip at the 99th quantile. na.rm because regions on contigs with no coverage are NA.
clip_val <- round(quantile(smooth_region_df$max_w_depth, probs=c(0.99), na.rm = TRUE))
smooth_region_df$max_w_depth[which(smooth_region_df$max_w_depth > clip_val)] <- clip_val
# and set anything in the 1% quantile to NA as uninformative
clip_val_low <- round(quantile(smooth_region_df$max_w_depth, probs=c(0.01), na.rm = TRUE))
smooth_region_df$max_w_depth[which(smooth_region_df$max_w_depth < clip_val_low)] <- NA

# ggplot(smooth_region_df, aes(x=start, y=max_w_depth)) +
#     geom_point() +
#     theme_bw() +
#     facet_wrap(.~chr)
#
# ggplot(smooth_region_df, aes(x=chr_i, y=max_w_depth)) +
#     geom_point() +
#     theme_bw() +
#     facet_wrap(.~chr)
#
# ggplot(smooth_region_df[smooth_region_df$chr=='chr1',], aes(x=start, y=max_w_depth)) +
#     geom_point() +
#     theme_bw()
# ggplot(smooth_region_df[smooth_region_df$chr=='chr1',], aes(x=chr_i, y=max_w_depth)) +
#     geom_point() +
#     theme_bw()
#

# write result
write.table(smooth_region_df, out_f, sep='\t', quote=F, row.names = F, col.names = T)
