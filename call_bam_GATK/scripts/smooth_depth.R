#!/usr/bin/env Rscript
# new LOH calling script on the project data
# working with the somatic variant calls, and depth from a single sample
library(ggplot2)
library(dplyr)

# args = commandArgs(trailingOnly=TRUE)
# if (length(args)!=3) {
#     stop("provide depth file, regions file, outfile as args 1, 2, 3")
# }

# depth_f <- gzfile(args[1], 'rt')
# regions_f <- args[2]
# out_f <- args[3]

depth_f <- gzfile(snakemake@input[[1]], 'rt')
regions_f <- snakemake@input[[2]]
out_f <- snakemake@output[[1]]

# read depth file
depth_df <- read.table(depth_f, header=F, sep='\t')
colnames(depth_df) <- c('chr', 'pos', 'depth')
depth_df$chr <- factor(depth_df$chr, levels=unique(depth_df$chr))

# read regions file
regions_df <- read.table(regions_f, sep='\t', quote='', header=F)
colnames(regions_df) <- c('chr', 'start', 'end', 'strand', 'name')
rownames(regions_df) <- regions_df$name
regions_df <- regions_df[order(as.numeric(gsub('chr', '', regions_df$chr)), regions_df$start),]
regions_df$chr <- factor(regions_df$chr, levels=unique(regions_df$chr))
regions_df$size <- regions_df$end - regions_df$start
regions_df$i <- 1:nrow(regions_df)
# within-chromosome index of each region. lapply, not sapply: when every chromosome has
# the same number of regions (e.g. exactly one, as with whole-contig targets) sapply
# simplifies to a vector and do.call(c, .) then fails with "second argument must be a list".
regions_df$chr_i <- unlist(lapply(table(regions_df$chr), function(x) seq_len(x)),
                           use.names = FALSE)

# smoothing function 
# sliding window with size w and step s
# have to deal with the fact that there's discontiguous blocks
smooth_depth_df <- function(depth_df, w=200, s=100){
    # first split into chr
    smooth_df_list <- lapply(unique(depth_df$chr), function(chr){
        print(chr)
        depth_df_chr <- depth_df[depth_df$chr==chr, ]
        keep_going <- TRUE
        smooth_df <- data.frame()
        start <- depth_df_chr[1, "pos"] 
        start_ind <- 1
        while(keep_going){
            w_stop <- start + w
            stop_ind <- start_ind + w
            if(stop_ind > nrow(depth_df_chr)){
                keep_going <- F
                break
            }
            stop_ind_pos <- depth_df_chr[stop_ind, "pos"]
            if (!is.na(stop_ind_pos) & (w_stop == stop_ind_pos)){
                mean_depth <- mean(depth_df_chr[start_ind:stop_ind, "depth"])
                new_df <- data.frame(chr=chr, start=start, end=w_stop, mean_depth=mean_depth)
                smooth_df <- rbind(smooth_df, new_df)    
            }
            start <- start+s 
            start_ind <- which(depth_df_chr[, "pos"] == start)
            if(length(start_ind)==0){
                # we've run past the end of a contiguous block
                # start at the next position with data
                start_ind <- which(depth_df_chr[, "pos"] > start)[1]
                start <- depth_df_chr[start_ind, "pos"]
            }
        }
        return(smooth_df)
        
    })
    return(smooth_df_list)
}


# could operate on each region, and find the max depth of a window of size w in that region 
smooth_depth_df_regions <- function(depth_df, regions_df, w=200, s=100){
    res_df <- regions_df
    res_df$max_w_depth <- 0
    
    depth_df$chr <- factor(depth_df$chr, levels=unique(depth_df$chr))
    dl <- group_split(depth_df, chr)
    names(dl) <- unique(depth_df$chr)
    current_r_chr <- ''
    for(i in 1:nrow(regions_df)){
    # for(i in 1:1000){
        r_chr <- regions_df[i, "chr"]
        if(current_r_chr != r_chr){
            print(as.character(r_chr))
            current_r_chr <- as.character(r_chr)
        }
        r_start <- regions_df[i, "start"]
        r_end <- regions_df[i, "end"]
        # A target contig can be missing from the depth table entirely: `samtools depth`
        # emits nothing for a contig with no alignments, even with -a. That is a real
        # situation (failed capture, a contig absent from the data), so record NA and
        # carry on rather than aborting the whole sample.
        if(!(as.character(r_chr) %in% names(dl))){
            res_df[i, 'max_w_depth'] <- NA
            next
        }
        depth_df_r <- dl[[as.character(r_chr)]]
        depth_df_r <- as.data.frame(depth_df_r[depth_df_r$pos >= r_start & depth_df_r$pos <= r_end, ])
        if(nrow(depth_df_r) == 0){
            res_df[i, 'max_w_depth'] <- NA
            next
        }
        # Regions shorter than the smoothing window cannot be windowed; summarise the
        # whole region instead of failing.
        if(nrow(depth_df_r) < w){
            res_df[i, 'max_w_depth'] <- round(mean(depth_df_r[, "depth"], na.rm = TRUE))
            next
        }
        max_w_depth <- 0
        start_ind <- 1
        end_ind <- start_ind + w
        while(end_ind <= nrow(depth_df_r)) {
            w_depth <- round(mean(depth_df_r[start_ind:end_ind, "depth"], na.rm = T))
            max_w_depth <- max(max_w_depth, w_depth)
            start_ind <- start_ind + s
            end_ind <- end_ind + s
        }
        res_df[i, 'max_w_depth'] <- max_w_depth
    }
    return(res_df)
}
    

smooth_region_df <- smooth_depth_df_regions(depth_df, regions_df, w=200, s=100)
smooth_region_df$chr <- factor(smooth_region_df$chr, levels=unique(smooth_region_df$chr))
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






