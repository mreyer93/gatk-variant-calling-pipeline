#!/usr/bin/env Rscript
# Regenerates the figures in example/README.md from a completed pipeline run.
# Run after test/run_test.sh:   Rscript example/make_figures.R test/data/results
suppressPackageStartupMessages({ library(ggplot2) })

res <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(res)) res <- "test/data/results"
figdir <- "example/figures"; outdir <- "example/results"
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

PAL <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860")
th <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey30", size = 10))
sv <- function(p, n, w = 8, h = 4.6) {
  ggsave(file.path(figdir, paste0(n, ".png")), p, width = w, height = h, dpi = 150, bg = "white")
  cat("  wrote", n, "\n")
}

# Read a Mutect2 VCF into a data frame, pulling AF and DP out of the FORMAT fields.
read_vcf <- function(path, label) {
  if (!file.exists(path)) return(NULL)
  ln <- readLines(path, warn = FALSE)
  hdr <- grep("^#CHROM", ln)
  body <- ln[(hdr + 1):length(ln)]
  body <- body[nzchar(body)]
  if (!length(body)) return(NULL)
  cols <- strsplit(sub("^#", "", ln[hdr]), "\t")[[1]]
  m <- do.call(rbind, strsplit(body, "\t"))
  d <- as.data.frame(m, stringsAsFactors = FALSE)
  colnames(d) <- cols
  fmt <- strsplit(d$FORMAT, ":")
  # last column is the sample of interest (tumour, for the paired calls)
  smp <- strsplit(d[[ncol(d)]], ":")
  getf <- function(k) mapply(function(f, s) {
      i <- match(k, f); if (is.na(i) || i > length(s)) NA_character_ else s[i] }, fmt, smp)
  d$AF <- suppressWarnings(as.numeric(sub(",.*$", "", getf("AF"))))
  d$DP <- suppressWarnings(as.numeric(getf("DP")))
  d$mode <- label
  d$pass <- d$FILTER == "PASS"
  d[, c("CHROM", "POS", "REF", "ALT", "FILTER", "AF", "DP", "mode", "pass")]
}

tvr_t <- read_vcf(file.path(res, "02_variants_reference/01_gatk_variant_calling/filtered/TUMOUR_filtered.vcf"), "Tumour vs reference")
tvr_n <- read_vcf(file.path(res, "02_variants_reference/01_gatk_variant_calling/filtered/NORMAL_filtered.vcf"), "Normal vs reference")
tvn   <- read_vcf(file.path(res, "03_variants_TvN/01_gatk_variant_calling/filtered/PATIENT_A_Dx_filtered.vcf"), "Tumour vs matched normal")
all_v <- do.call(rbind, Filter(Negate(is.null), list(tvr_n, tvr_t, tvn)))
write.table(all_v, file.path(outdir, "variant_calls.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- 1. calls by mode ---------------------------------------------------------
cnt <- as.data.frame(table(all_v$mode, ifelse(all_v$pass, "PASS", "filtered out")))
colnames(cnt) <- c("mode", "status", "n")
cnt$mode <- factor(cnt$mode, levels = c("Normal vs reference", "Tumour vs reference",
                                        "Tumour vs matched normal"))
sv(ggplot(cnt, aes(mode, n, fill = status)) +
     geom_col() +
     geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.4, colour = "white") +
     scale_fill_manual(values = c(PASS = PAL[3], `filtered out` = "grey65"), name = NULL) +
     labs(title = "Subtracting the matched normal removes germline variants",
          subtitle = "same tumour sample, called three ways",
          x = NULL, y = "Variant calls") + th, "calls_by_mode", 8, 4.4)

# --- 2. VAF distribution -------------------------------------------------------
v <- all_v[is.finite(all_v$AF) & all_v$pass, ]
if (nrow(v)) {
  sv(ggplot(v, aes(AF, fill = mode)) +
       geom_histogram(bins = 25, colour = NA, alpha = 0.85) +
       facet_wrap(~mode, ncol = 1, scales = "free_y") +
       scale_fill_manual(values = PAL, guide = "none") +
       labs(title = "Variant allele fraction of PASS calls",
            subtitle = "germline variants cluster near 0.5 and 1.0; somatic calls sit lower",
            x = "Allele fraction", y = "Variants") + th, "vaf_distribution", 8, 6)
}

# --- 3. filter reasons ---------------------------------------------------------
fr <- all_v[!all_v$pass, ]
if (nrow(fr)) {
  reasons <- unlist(strsplit(fr$FILTER, ";"))
  rd <- as.data.frame(table(reasons)); colnames(rd) <- c("reason", "n")
  rd <- rd[order(-rd$n), ]
  sv(ggplot(head(rd, 10), aes(reorder(reason, n), n)) +
       geom_col(fill = PAL[4]) + coord_flip() +
       labs(title = "Why calls were filtered out",
            subtitle = "GATK FilterMutectCalls annotations, all three calling modes",
            x = NULL, y = "Calls") + th, "filter_reasons", 8, 4)
  write.table(rd, file.path(outdir, "filter_reasons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# --- 4. depth ------------------------------------------------------------------
dp <- file.path(res, "01_prepare_bam/depth_smoothed_join.tsv")
if (file.exists(dp)) {
  d <- read.delim(dp, check.names = FALSE)
  meta <- intersect(c("chr","start","end","strand","name","size","i","chr_i"), colnames(d))
  samples <- setdiff(colnames(d), meta)
  dd <- data.frame(sample = samples,
                   depth = as.numeric(d[1, samples]), stringsAsFactors = FALSE)
  write.table(dd, file.path(outdir, "depth.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  sv(ggplot(dd, aes(sample, depth, fill = sample)) + geom_col(width = 0.6) +
       geom_text(aes(label = round(depth)), vjust = -0.5, size = 3.6) +
       scale_fill_manual(values = PAL, guide = "none") +
       expand_limits(y = max(dd$depth) * 1.15) +
       labs(title = "On-target sequencing depth", subtitle = "smoothed max depth over the target region",
            x = NULL, y = "Depth (x)") + th, "depth", 6, 4)
}
cat("done\n")
