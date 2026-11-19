# Copy number of the SEQC2 HCC1395 exome: the pipeline's FACETS run against SEQC2's own
# copy-number calls on the same data.
#
# References (downloaded by benchmark/seqc2/download_seqc2.sh into <seqc2 dir>/cnv/):
#   - SEQC2's FACETS segment tables for the 12 exome replicates (WES_*_fit_cncf.txt); WES_IL_1
#     is the same Illumina-centre library pair the pipeline run uses
#   - SEQC2's ASCAT segments for the same tumour library (WES_IL_T_1.copynumber.caveman.csv:
#     segment, chromosome, start, end, normal total and minor copy number, tumour total and
#     minor copy number)
#   - the karyotype and copy-number findings in Fang et al. 2021 (Nat Biotechnol 39:1151,
#     https://doi.org/10.1038/s41587-021-00993-6): tumour purity 99% from copy-number
#     analysis, a hypertetraploid karyotype (64 to 79 chromosomes plus 38 to 63 marker
#     chromosomes; Extended Data Fig. 5a), clonal copy-neutral LOH of 17q, and loss of one
#     copy of 6p, 16q and chrX in the normal line HCC1395BL (Extended Data Figs. 5b and 6)
#
# Run from the repository root in the pipeline environment, after the somatic pipeline has
# run with run_facets: True:
#   Rscript benchmark/seqc2/compare_facets.R <run output_directory> <seqc2 dir> <results dir>

options(width = 120)
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

# Step 1. Inputs
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript compare_facets.R <run output_directory> <seqc2 dir> <results dir>")
run.dir <- normalizePath(args[1])
cnv.dir <- file.path(normalizePath(args[2]), "cnv")
out.dir <- args[3]
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
out.dir <- normalizePath(out.dir)
purity.f <- file.path(run.dir, "04_FACETS/HCC1395_Dx_purity.txt")
cncf.f <- file.path(run.dir, "04_FACETS/HCC1395_Dx_cncf.tsv")
ascat.f <- file.path(cnv.dir, "ascat/WES_IL_T_1.copynumber.caveman.csv")
seqc2.fs <- list.files(file.path(cnv.dir, "facets"), pattern = "_fit_cncf.txt$", full.names = TRUE)
for (f in c(purity.f, cncf.f, ascat.f)) if (!file.exists(f)) stop("Missing input: ", f)
if (length(seqc2.fs) == 0) stop("No SEQC2 FACETS tables in ", file.path(cnv.dir, "facets"))

# Segment tables in one shape: chromosome (1-22, 23 = X), start, end, total and minor copy
# number, cellular fraction, and for FACETS the clustered log ratio and SNP count
chrom.num <- function(x) { x <- sub("^chr", "", x); x[x == "X"] <- "23"; as.integer(x) }
pipe <- read.delim(cncf.f)
fit <- read.delim(purity.f)
facets.segs <- function(d) with(d, data.frame(chrom = chrom.num(chrom), start, end, tcn = tcn.em, lcn = lcn.em,
                                              cf = cf.em, cnlr = cnlr.median.clust, marks = num.mark))
segs.pipe <- facets.segs(pipe)
seqc2 <- lapply(seqc2.fs, function(f) facets.segs(read.delim(f)))   # first column of these files is a row index
names(seqc2) <- sub("_fit_cncf.txt$", "", basename(seqc2.fs))
a <- read.csv(ascat.f, header = FALSE, col.names = c("id", "chr", "start", "end", "n_tcn", "n_lcn", "tcn", "lcn"))
segs.ascat <- data.frame(chrom = chrom.num(a$chr), start = a$start, end = a$end, tcn = a$tcn, lcn = a$lcn,
                         cf = NA, cnlr = NA, marks = NA)

# Step 2. Purity and ploidy
# FACETS reports purity and ploidy with its fit; SEQC2 published only the segment tables, so
# both are recovered from them, the same way for every FACETS table:
#   purity: the most common cellular fraction (by length) among segments that are not 2 copies
#     with one of each allele, since clonal changes carry the purity as their cellular fraction
#   ploidy: FACETS' own definition, (2^(1 - dipLogR) - 2(1 - purity)) / purity, where dipLogR
#     (the log ratio of the two-copy state) is the SNP-weighted median over autosomal segments of
#     the clustered log ratio minus log2(1 - cf + cf * tcn / 2), FACETS' expected shift
# The recovery is checked against the pipeline's reported purity, dipLogR and ploidy. Mean total
# copy number (length-weighted) is also given, as the only ploidy measure ASCAT's table allows.
wmed <- function(x, w) { o <- order(x); x <- x[o]; w <- w[o]; x[which(cumsum(w) >= sum(w) / 2)[1]] }
derive <- function(s) {
    s$len <- as.numeric(s$end - s$start)
    ok <- !is.na(s$tcn)
    changed <- ok & !is.na(s$cf) & !(s$tcn == 2 & s$lcn %in% 1)
    purity <- if (any(changed)) as.numeric(names(which.max(tapply(s$len[changed], round(s$cf[changed], 3), sum)))) else NA
    k <- ok & !is.na(s$cf) & !is.na(s$cnlr) & s$tcn > 0 & s$chrom <= 22
    dip <- if (any(k)) wmed(s$cnlr[k] - log2(1 - s$cf[k] + s$cf[k] * s$tcn[k] / 2), s$marks[k]) else NA
    data.frame(purity = purity, dipLogR = dip, ploidy = (2^(1 - dip) - 2 * (1 - purity)) / purity,
               mean_total_cn = sum(s$tcn[ok] * s$len[ok]) / sum(s$len[ok]),
               loh_fraction = sum(s$len[ok & s$lcn %in% 0]) / sum(s$len[ok & !is.na(s$lcn)]),
               segments = nrow(s))
}
pp <- bind_rows(
    cbind(source = "pipeline FACETS, as reported", data.frame(purity = fit$purity, dipLogR = fit$dipLogR, ploidy = fit$ploidy)),
    cbind(source = "pipeline FACETS, recovered", derive(segs.pipe)),
    cbind(source = "SEQC2 FACETS WES_IL_1 (same pair)", derive(seqc2[["WES_IL_1"]])),
    cbind(source = "SEQC2 ASCAT WES_IL_T_1 (same tumour)", derive(segs.ascat)))
cat("\nStep 2. Purity, ploidy and LOH fraction\n")
print(pp, row.names = FALSE, digits = 3)
reps <- bind_rows(lapply(names(seqc2), function(r) cbind(replicate = r, derive(seqc2[[r]]))))
cat("\n  SEQC2 FACETS, all exome replicates\n")
print(reps, row.names = FALSE, digits = 3)
write.table(bind_rows(pp, reps %>% rename(source = replicate) %>% mutate(source = paste("SEQC2 FACETS", source))),
            file.path(out.dir, "facets_purity_ploidy.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

####
# Step 2. Purity, ploidy and LOH fraction
#                                source purity dipLogR ploidy mean_total_cn loh_fraction segments
#          pipeline FACETS, as reported  0.948  -0.611   3.11            NA           NA       NA
#            pipeline FACETS, recovered  0.952  -0.610   3.11          3.32        0.585      251
#     SEQC2 FACETS WES_IL_1 (same pair)  0.937  -0.635   3.18          3.36        0.595      341
#  SEQC2 ASCAT WES_IL_T_1 (same tumour)     NA      NA     NA          3.84        0.659     6253
#
#   SEQC2 FACETS, all exome replicates
#  replicate purity dipLogR ploidy mean_total_cn loh_fraction segments
#   WES_EA_1  0.848  -1.470   6.18          5.91        0.592      387
#   WES_FD_1  0.831  -1.391   5.91          5.85        0.587      255
#   WES_FD_2  0.911  -0.575   3.08          3.26        0.585      255
#   WES_FD_3  0.889  -0.575   3.10          3.28        0.606      269
#   WES_IL_1  0.937  -0.635   3.18          3.36        0.595      341
#   WES_IL_2  0.871  -0.640   3.28          3.47        0.596      689
#   WES_IL_3  0.948  -0.591   3.07          3.21        0.611      501
#   WES_LL_1  0.902  -0.579   3.10          3.19        0.603      266
#   WES_NC_1  0.896  -1.345   5.44          5.44        0.593      348
#   WES_NV_1  0.944  -0.613   3.12          3.33        0.598      659
#   WES_NV_2  0.952  -0.595   3.07          3.27        0.611      712
#   WES_NV_3  0.958  -0.600   3.08          3.28        0.608      699
####

# Step 3. Segment-by-segment agreement
# Share of the genome (bases covered by both segmentations, autosomes) where two calls agree on
# total copy number, on minor copy number, and on LOH (minor copy number 0). Total copy
# number only agrees when the two fits chose the same ploidy, so LOH is the comparison that
# survives a different ploidy solution. 6p and 16q are reported separately because the normal
# line has lost a copy of each (below).
normal.loss <- function(s) (s$chrom == 6 & s$end <= 60e6) | (s$chrom == 16 & s$start >= 36e6)
agree <- function(x, y) {
    # interval intersection of two segment tables, one row per overlapping piece
    out <- list()
    for (ch in intersect(unique(x$chrom), unique(y$chrom))) {
        xx <- x[x$chrom == ch, ]; yy <- y[y$chrom == ch, ]
        for (i in seq_len(nrow(xx))) {
            j <- which(yy$start < xx$end[i] & yy$end > xx$start[i])
            if (!length(j)) next
            s <- pmax(xx$start[i], yy$start[j]); e <- pmin(xx$end[i], yy$end[j])
            out[[length(out) + 1]] <- data.frame(chrom = ch, start = s, end = e,
                tcn.x = xx$tcn[i], lcn.x = xx$lcn[i], tcn.y = yy$tcn[j], lcn.y = yy$lcn[j])
        }
    }
    o <- bind_rows(out)
    o$len <- as.numeric(o$end - o$start)
    o
}
summarise.agree <- function(o, label) {
    o <- o[o$chrom <= 22 & !is.na(o$tcn.x) & !is.na(o$tcn.y), ]
    bind_rows(lapply(c("autosomes", "autosomes outside 6p and 16q"), function(scope) {
        k <- if (scope == "autosomes") rep(TRUE, nrow(o)) else !normal.loss(o)
        m <- o[k & !is.na(o$lcn.x) & !is.na(o$lcn.y), ]
        data.frame(comparison = label, scope = scope, shared_Mb = round(sum(o$len[k]) / 1e6, 1),
                   same_total_cn = sum(o$len[k][o$tcn.x[k] == o$tcn.y[k]]) / sum(o$len[k]),
                   same_minor_cn = sum(m$len[m$lcn.x == m$lcn.y]) / sum(m$len),
                   same_loh_call = sum(m$len[(m$lcn.x == 0) == (m$lcn.y == 0)]) / sum(m$len))
    }))
}
ov.seqc2 <- agree(segs.pipe, seqc2[["WES_IL_1"]])
ov.ascat <- agree(segs.pipe, segs.ascat)
ov.ref <- agree(seqc2[["WES_IL_1"]], segs.ascat)
conc <- bind_rows(summarise.agree(ov.seqc2, "pipeline vs SEQC2 FACETS WES_IL_1"),
                  summarise.agree(ov.ascat, "pipeline vs SEQC2 ASCAT WES_IL_T_1"),
                  summarise.agree(ov.ref, "SEQC2 FACETS vs SEQC2 ASCAT (reference spread)"))
cat("\nStep 3. Agreement, as a share of bases covered by both\n")
print(conc %>% mutate(across(starts_with("same"), ~ round(.x, 3))), row.names = FALSE)
write.table(conc, file.path(out.dir, "facets_agreement.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

####
# Step 3. Agreement, as a share of bases covered by both
#                                      comparison                        scope shared_Mb same_total_cn same_minor_cn
#               pipeline vs SEQC2 FACETS WES_IL_1                    autosomes    2746.3         0.870         0.945
#               pipeline vs SEQC2 FACETS WES_IL_1 autosomes outside 6p and 16q    2656.6         0.880         0.945
#              pipeline vs SEQC2 ASCAT WES_IL_T_1                    autosomes    2709.9         0.296         0.807
#              pipeline vs SEQC2 ASCAT WES_IL_T_1 autosomes outside 6p and 16q    2612.8         0.285         0.807
#  SEQC2 FACETS vs SEQC2 ASCAT (reference spread)                    autosomes    2697.9         0.290         0.803
#  SEQC2 FACETS vs SEQC2 ASCAT (reference spread) autosomes outside 6p and 16q    2601.1         0.283         0.803
#  same_loh_call
#          0.999
#          0.999
#          0.917
#          0.917
#          0.925
#          0.926
####

# Step 4. Regions named by SEQC2
# 17q: clonal copy-neutral LOH in the tumour (Fang et al. 2021). 6p, 16q and X: one copy lost
# in the normal line, so FACETS (which takes the normal as two copies) sees relatively more
# tumour DNA there and has no heterozygous normal SNPs to read allele balance from. Arm
# boundaries are approximate (hg38 centromeres lie near 58-62 Mb on chr6, 36-38 Mb on chr16 and
# 23-27 Mb on chr17); the states reported are the three largest by length in each region.
region.cn <- function(s, ch, from, to) {
    s <- s[s$chrom == ch & s$end > from & s$start < to & !is.na(s$tcn), ]
    s$len <- pmin(s$end, to) - pmax(s$start, from)
    s %>% group_by(tcn, lcn) %>% summarise(Mb = round(sum(len) / 1e6, 1), .groups = "drop") %>%
        arrange(desc(Mb)) %>% mutate(state = paste0(tcn, "/", ifelse(is.na(lcn), "NA", lcn), " (", Mb, " Mb)")) %>%
        pull(state) %>% head(3) %>% paste(collapse = "; ")
}
named <- data.frame(region = c("17q (tumour cn-LOH)", "6p (normal lost a copy)", "16q (normal lost a copy)", "X (normal lost a copy)"),
                    chrom = c(17, 6, 16, 23), from = c(27.4e6, 0, 36.8e6, 0), to = c(84e6, 59.8e6, 91e6, 157e6))
named$pipeline <- mapply(region.cn, list(segs.pipe), named$chrom, named$from, named$to)
named$seqc2_facets <- mapply(region.cn, list(seqc2[["WES_IL_1"]]), named$chrom, named$from, named$to)
named$ascat <- mapply(region.cn, list(segs.ascat), named$chrom, named$from, named$to)
cat("\nStep 4. Largest total/minor copy-number states in the named regions\n")
for (i in seq_len(nrow(named))) {
    cat(sprintf("  %s\n    pipeline FACETS: %s\n    SEQC2 FACETS:    %s\n    SEQC2 ASCAT:     %s\n",
                named$region[i], named$pipeline[i], named$seqc2_facets[i], named$ascat[i]))
}

# Figure 3. Genome-wide total and minor copy number from the three calls (numbered after
# evaluate_seqc2.R's two figures, as in manual/benchmark.md)
lens <- bind_rows(segs.pipe, seqc2[["WES_IL_1"]], segs.ascat) %>% filter(chrom <= 23) %>%
    group_by(chrom) %>% summarise(len = max(end), .groups = "drop") %>% arrange(chrom) %>%
    mutate(offset = cumsum(as.numeric(lag(len, default = 0))), mid = offset + len / 2)
prep <- function(s, label) {
    s %>% filter(chrom <= 23, !is.na(tcn)) %>% left_join(lens, by = "chrom") %>%
        mutate(x0 = start + offset, x1 = end + offset, call = label)
}
fd <- bind_rows(prep(segs.pipe, sprintf("Pipeline FACETS (purity %.2f, ploidy %.2f)", fit$purity, fit$ploidy)),
                prep(seqc2[["WES_IL_1"]], "SEQC2 FACETS, same library pair"),
                prep(segs.ascat, "SEQC2 ASCAT, same tumour library"))
fd$call <- factor(fd$call, levels = unique(fd$call))
cap <- 8   # copy numbers above 8 drawn at 8
cn.pal <- c("total copy number" = "#2a78d6", "minor copy number" = "#eb6834")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid.col <- "#e6e5e1"; surface <- "#fcfcfb"
fig1 <- ggplot(fd) +
    geom_vline(xintercept = lens$offset[-1], colour = grid.col, linewidth = 0.3) +
    geom_segment(aes(x = x0, xend = x1, y = pmin(tcn, cap), yend = pmin(tcn, cap), colour = "total copy number"),
                 linewidth = 1.2, lineend = "butt") +
    geom_segment(data = fd[!is.na(fd$lcn), ], aes(x = x0, xend = x1, y = pmin(lcn, cap) - 0.12,
                 yend = pmin(lcn, cap) - 0.12, colour = "minor copy number"), linewidth = 1.2, lineend = "butt") +
    facet_wrap(~call, ncol = 1) +
    scale_colour_manual(values = cn.pal, breaks = names(cn.pal)) +
    scale_x_continuous(breaks = lens$mid, labels = c(1:22, "X"), expand = c(0.005, 0)) +
    scale_y_continuous(breaks = 0:cap, limits = c(-0.3, cap + 0.3)) +
    labs(x = NULL, y = "Copy number", title = "HCC1395 copy number, three calls on the same exome data",
         subtitle = "Minor copy number drawn slightly below total where they coincide; values above 8 drawn at 8") +
    theme_minimal(base_size = 10) + theme(
        plot.background = element_rect(fill = surface, colour = NA),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = grid.col, linewidth = 0.3),
        axis.text = element_text(colour = ink2), axis.text.x = element_text(size = 7),
        axis.title = element_text(colour = ink2), strip.text = element_text(colour = ink, hjust = 0),
        plot.title = element_text(colour = ink, face = "bold"), plot.subtitle = element_text(colour = ink2),
        legend.position = "top", legend.justification = "left", legend.title = element_blank(),
        legend.text = element_text(colour = ink2), plot.title.position = "plot")
ggsave(file.path(out.dir, "figure3_copy_number.png"), fig1, width = 10, height = 6.5, dpi = 150, bg = surface)

####
# Step 4. Largest total/minor copy-number states in the named regions
#   17q (tumour cn-LOH)
#     pipeline FACETS: 2/0 (39.8 Mb); 5/0 (8.3 Mb); 4/0 (7.1 Mb)
#     SEQC2 FACETS:    2/0 (39.6 Mb); 6/0 (8.3 Mb); 4/0 (7.1 Mb)
#     SEQC2 ASCAT:     3/0 (33.4 Mb); 5/0 (10.3 Mb); 4/0 (5.3 Mb)
#   6p (normal lost a copy)
#     pipeline FACETS: 7/NA (30.8 Mb); 4/NA (13.2 Mb); 9/NA (9.2 Mb)
#     SEQC2 FACETS:    8/NA (22.7 Mb); 7/NA (14.8 Mb); 5/NA (13.1 Mb)
#     SEQC2 ASCAT:     7/0 (40.9 Mb); 4/0 (8.8 Mb); 8/2 (5.4 Mb)
#   16q (normal lost a copy)
#     pipeline FACETS: 6/NA (23.4 Mb); 5/NA (19 Mb); 8/NA (11 Mb)
#     SEQC2 FACETS:    6/NA (24.5 Mb); 9/NA (11 Mb); 5/NA (8.2 Mb)
#     SEQC2 ASCAT:     6/1 (34.8 Mb); 7/1 (8.1 Mb); 8/0 (0.9 Mb)
#   X (normal lost a copy)
#     pipeline FACETS: 2/NA (77.9 Mb); 1/NA (68.4 Mb); 4/NA (3.1 Mb)
#     SEQC2 FACETS:    2/NA (78.1 Mb); 1/NA (67.3 Mb); 4/NA (3.2 Mb)
#     SEQC2 ASCAT:     5/0 (76.7 Mb); 3/0 (24 Mb); 2/0 (19.2 Mb)
####
