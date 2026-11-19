# Somatic accuracy of the tumour-vs-normal calls on the SEQC2 HCC1395 whole exome.
#
# Scores a finished run of call_bam_GATK (tumour HCC1395 vs normal HCC1395BL, Illumina-centre
# library WES_IL_T_1 / WES_IL_N_1) against the SEQC2 somatic reference call set v1.2.1
# (HighConf + MedConf SNVs and indels; Fang et al. 2021, Nat Biotechnol 39:1151,
# https://doi.org/10.1038/s41587-021-00993-6). Scoring follows the SEQC2 WES benchmark
# (Xiao et al. 2021, Nat Biotechnol 39:1141, https://doi.org/10.1038/s41587-021-00994-5):
# calls and truth restricted to the exome capture targets, recall = TP/(TP+FN),
# precision = TP/(TP+FP).
#
# Matching: both sides are normalised the same way (multiallelic records split, indels
# left-aligned, Mutect2's multi-base substitutions split into single SNVs), then matched on
# chromosome, position, REF and ALT. A record is in a region when its POS is.
#
# Run from the repository root in the pipeline environment (envs/environment.yml):
#   Rscript benchmark/seqc2/evaluate_seqc2.R <run output_directory> <reference.fasta> \
#       <seqc2 data dir> <results dir> [calls.vcf]
# The seqc2 data dir is laid out by benchmark/seqc2/download_seqc2.sh. The run must use
# benchmark/seqc2/bam_metadata.tsv (samples TUMOUR and NORMAL, patient HCC1395, timepoint
# Dx). The optional calls.vcf scores a re-filtered version of the run's calls instead of
# its own FilterMutectCalls output (benchmark/seqc2/backtest_filters.sh uses this).

options(width = 120)
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

# Step 1. Inputs
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 4:5) {
    stop("Usage: Rscript evaluate_seqc2.R <run output_directory> <reference.fasta> <seqc2 dir> <results dir> [calls.vcf]")
}
run.dir <- normalizePath(args[1])
ref.fa <- normalizePath(args[2])
seqc2.dir <- normalizePath(args[3])
out.dir <- args[4]
dir.create(file.path(out.dir, "work"), recursive = TRUE, showWarnings = FALSE)
out.dir <- normalizePath(out.dir)
w <- function(f) file.path(out.dir, "work", f)

calls.vcf <- if (length(args) == 5) normalizePath(args[5]) else
    file.path(run.dir, "03_variants_TvN/01_gatk_variant_calling/filtered/HCC1395_Dx_filtered.vcf")
tumour.bam <- file.path(run.dir, "01_prepare_bam/recalibrate/TUMOUR.fixchr.bam")
normal.bam <- file.path(run.dir, "01_prepare_bam/recalibrate/NORMAL.fixchr.bam")
targets.bed <- file.path(seqc2.dir, "S07604624_Covered_human_all_v6_plus_UTR.liftover.to.hg38.bed6.gz")
hc.bed <- file.path(seqc2.dir, "truth/High-Confidence_Regions_v1.2.bed")
truth.snv <- file.path(seqc2.dir, "truth/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz")
truth.indel <- file.path(seqc2.dir, "truth/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz")
superset.snv <- file.path(seqc2.dir, "truth/sSNV.MSDUKT.superSet.v1.2.vcf.gz")
superset.indel <- file.path(seqc2.dir, "truth/sINDEL.MDKT.superSet.v1.2.vcf.gz")
for (f in c(calls.vcf, tumour.bam, normal.bam, targets.bed, hc.bed, truth.snv, truth.indel,
            superset.snv, superset.indel, ref.fa)) {
    if (!file.exists(f)) stop("Missing input: ", f)
}

# shell commands run through bash with pipefail, so a failure anywhere stops the script
sh <- function(cmd) {
    status <- system2("bash", c("-o", "pipefail", "-c", shQuote(cmd)))
    if (status != 0) stop("Command failed (", status, "): ", cmd)
}
read.tsv <- function(f, cols) {
    if (file.size(f) == 0) {
        return(setNames(data.frame(matrix(nrow = 0, ncol = length(cols))), cols))
    }
    read.delim(f, header = FALSE, col.names = cols, colClasses = "character", quote = "",
               na.strings = ".")
}
bed.mb <- function(f) {
    b <- read.delim(f, header = FALSE)
    sum(as.numeric(b$V3 - b$V2)) / 1e6
}

# Step 2. Evaluation regions
# Capture targets (primary chromosomes), the padded region the pipeline calls in (targets
# +/- 100 bp, as in V6_UTR_w100.bed), and the SEQC2 high-confidence regions. The scored
# region is targets n high-confidence; the padded version is a secondary check.
primary <- "^chr([0-9]+|X|Y)\t"
sh(sprintf("gzip -dc %s | grep -E '%s' | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > %s",
           shQuote(targets.bed), primary, w("targets.bed")))
sh(sprintf("awk 'BEGIN{OFS=\"\\t\"} {s=$2-100; if (s<0) s=0; print $1, s, $3+100}' %s | sort -k1,1 -k2,2n | bedtools merge > %s",
           w("targets.bed"), w("targets_w100.bed")))
sh(sprintf("grep -E '%s' %s | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge > %s",
           primary, shQuote(hc.bed), w("hc.bed")))
sh(sprintf("bedtools intersect -a %s -b %s | sort -k1,1 -k2,2n | bedtools merge > %s", w("targets.bed"), w("hc.bed"), w("eval.bed")))
sh(sprintf("bedtools intersect -a %s -b %s | sort -k1,1 -k2,2n | bedtools merge > %s", w("targets_w100.bed"), w("hc.bed"), w("eval_w100.bed")))
regions <- data.frame(region = c("targets", "targets_w100", "hc", "eval", "eval_w100"))
regions$Mb <- sapply(regions$region, function(r) bed.mb(w(paste0(r, ".bed"))))
cat("\nStep 2. Region sizes (Mb)\n")
print(regions, row.names = FALSE, digits = 4)

####
# Step 2. Region sizes (Mb)
#        region      Mb
#       targets   90.54
#  targets_w100  133.99
#            hc 2485.47
#          eval   80.71
#     eval_w100  119.09
####

# Step 3. Normalise truth and calls
# Same normalisation for both: split multiallelics, left-align against the reference,
# split multi-base substitutions into SNVs; then drop exact duplicates.
norm <- sprintf("bcftools norm -f %s -m -any --atomize --old-rec-tag ORIGINAL --check-ref w -Ou 2>/dev/null | bcftools norm -d exact -Ou 2>/dev/null",
                shQuote(ref.fa))
for (t in c("snv", "indel")) {
    vcf <- if (t == "snv") truth.snv else truth.indel
    sh(sprintf("bcftools view %s -Ou | %s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%FILTER\\t%%INFO/TVAF\\n' > %s",
               shQuote(vcf), norm, w(paste0("truth_", t, ".tsv"))))
    sh(sprintf("bcftools view %s -Ou | %s | bcftools view -T %s -Ou | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' > %s",
               shQuote(vcf), norm, w("eval.bed"), w(paste0("truth_", t, "_in_eval.tsv"))))
    sh(sprintf("bcftools view %s -Ou | %s | bcftools view -T %s -Ou | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\n' > %s",
               shQuote(vcf), norm, w("eval_w100.bed"), w(paste0("truth_", t, "_in_eval_w100.tsv"))))
}
key <- function(d) paste(d$CHROM, d$POS, d$REF, d$ALT, sep = ":")
truth <- bind_rows(
    read.tsv(w("truth_snv.tsv"), c("CHROM", "POS", "REF", "ALT", "TRUTH_FILTER", "TVAF")) %>% mutate(type = "SNV"),
    read.tsv(w("truth_indel.tsv"), c("CHROM", "POS", "REF", "ALT", "TRUTH_FILTER", "TVAF")) %>% mutate(type = "indel"))
truth$key <- key(truth)
truth$TVAF <- as.numeric(truth$TVAF)
truth$in_eval <- truth$key %in% c(key(read.tsv(w("truth_snv_in_eval.tsv"), c("CHROM", "POS", "REF", "ALT"))),
                                  key(read.tsv(w("truth_indel_in_eval.tsv"), c("CHROM", "POS", "REF", "ALT"))))
truth$in_eval_w100 <- truth$key %in% c(key(read.tsv(w("truth_snv_in_eval_w100.tsv"), c("CHROM", "POS", "REF", "ALT"))),
                                       key(read.tsv(w("truth_indel_in_eval_w100.tsv"), c("CHROM", "POS", "REF", "ALT"))))

# Every record FilterMutectCalls wrote, whatever its FILTER, so missed truth calls can be
# traced to the filter that removed them. Sample columns are found from the header.
hdr <- system2("bcftools", c("view", "-h", shQuote(calls.vcf)), stdout = TRUE)
tumour.name <- sub("^##tumor_sample=", "", grep("^##tumor_sample=", hdr, value = TRUE))
normal.name <- sub("^##normal_sample=", "", grep("^##normal_sample=", hdr, value = TRUE))
fmt <- "%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO/TLOD\\t%INFO/ORIGINAL\\t%INFO/POPAF[\\t%AD\\t%AF\\t%DP\\t%F1R2\\t%F2R1]\\n"
call.cols <- c("CHROM", "POS", "REF", "ALT", "FILTER", "TLOD", "ORIGINAL", "POPAF",
               "T_AD", "T_AF", "T_DP", "T_F1R2", "T_F2R1", "N_AD", "N_AF", "N_DP", "N_F1R2", "N_F2R1")
for (r in c("all", "eval", "eval_w100", "targets", "targets_w100")) {
    sel <- if (r == "all") "" else sprintf("| bcftools view -T %s -Ou", w(paste0(r, ".bed")))
    sh(sprintf("bcftools view -s %s,%s %s -Ou | %s %s | bcftools query -f '%s' > %s",
               tumour.name, normal.name, shQuote(calls.vcf), norm, sel, fmt, w(paste0("calls_", r, ".tsv"))))
}
calls <- read.tsv(w("calls_all.tsv"), call.cols)
calls$key <- key(calls)
for (r in c("eval", "eval_w100", "targets", "targets_w100")) {
    calls[[paste0("in_", r)]] <- calls$key %in% key(read.tsv(w(paste0("calls_", r, ".tsv")), call.cols))
}
calls$type <- ifelse(nchar(calls$REF) == 1 & nchar(calls$ALT) == 1, "SNV", "indel")
alt.of <- function(x) as.numeric(sub("^[^,]*,", "", x))   # AD is "ref,alt" after splitting
calls <- calls %>% mutate(T_ALT = alt.of(T_AD), N_ALT = alt.of(N_AD), T_AF = as.numeric(T_AF),
                          N_AF = as.numeric(N_AF), T_DP = as.numeric(T_DP), N_DP = as.numeric(N_DP),
                          TLOD = as.numeric(TLOD), POPAF = as.numeric(POPAF))
cat("\nStep 3. Normalised records\n")
cat(sprintf("  truth: %d SNVs and %d indels genome-wide; %d SNVs and %d indels in the scored region\n",
            sum(truth$type == "SNV"), sum(truth$type == "indel"),
            sum(truth$type == "SNV" & truth$in_eval), sum(truth$type == "indel" & truth$in_eval)))
cat(sprintf("  calls: %d records from FilterMutectCalls (%d PASS); %d PASS in the scored region\n",
            nrow(calls), sum(calls$FILTER == "PASS"), sum(calls$FILTER == "PASS" & calls$in_eval)))

####
# Step 3. Normalised records
#   truth: 39447 SNVs and 1625 indels genome-wide; 1159 SNVs and 48 indels in the scored region
#   calls: 14765 records from FilterMutectCalls (3648 PASS); 1108 PASS in the scored region
####

# Step 4. Precision and recall
# Two call sets: FilterMutectCalls PASS, and the pipeline's own keep rule for the annotated
# tables (FILTER exactly PASS or normal_artifact, AF >= 0.02 and DP >= 5 in either sample;
# call_bam_GATK/scripts/combine_annotation_TvN.R, without its CADD cut, which ranks
# functional impact rather than calling). 95% intervals are Wilson score intervals.
wilson <- function(x, n) {
    if (n == 0) return(c(NA, NA))
    p <- x / n; z <- 1.96
    centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
    c(centre - half, centre + half)
}
calls$keep_pass <- calls$FILTER == "PASS"
calls$keep_pipeline <- calls$FILTER %in% c("PASS", "normal_artifact") &
    (calls$T_AF >= 0.02 | calls$N_AF >= 0.02) & (calls$T_DP >= 5 | calls$N_DP >= 5)
score <- function(callset, region, vtype) {
    tr <- truth$key[truth$type == vtype & truth[[paste0("in_", region)]]]
    cl <- calls$key[calls$type == vtype & calls[[callset]] & calls[[paste0("in_", region)]]]
    tp <- sum(cl %in% tr); fp <- sum(!cl %in% tr); fn <- sum(!tr %in% cl)
    rc <- wilson(tp, tp + fn); pr <- wilson(tp, tp + fp)
    data.frame(callset = sub("keep_", "", callset), region = region, type = vtype, truth = length(tr),
               calls = length(cl), TP = tp, FP = fp, FN = fn,
               precision = tp / (tp + fp), prec_lo = pr[1], prec_hi = pr[2],
               recall = tp / (tp + fn), rec_lo = rc[1], rec_hi = rc[2],
               F1 = 2 * tp / (2 * tp + fp + fn))
}
scores <- bind_rows(lapply(c("keep_pass", "keep_pipeline"), function(cs)
    bind_rows(lapply(c("eval", "eval_w100"), function(r)
        bind_rows(lapply(c("SNV", "indel"), function(t) score(cs, r, t)))))))
cat("\nStep 4. Precision and recall (eval = targets n HC; eval_w100 = padded targets n HC)\n")
print(scores %>% mutate(across(precision:F1, ~ round(.x, 3))), row.names = FALSE)
write.table(scores, file.path(out.dir, "scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

####
# Step 4. Precision and recall (eval = targets n HC; eval_w100 = padded targets n HC)
#   callset    region  type truth calls   TP FP  FN precision prec_lo prec_hi recall rec_lo rec_hi    F1
#      pass      eval   SNV  1159  1056 1047  9 112     0.991   0.984   0.996  0.903  0.885  0.919 0.945
#      pass      eval indel    48    52   41 11   7     0.788   0.660   0.878  0.854  0.728  0.928 0.820
#      pass eval_w100   SNV  1666  1427 1402 25 264     0.982   0.974   0.988  0.842  0.823  0.858 0.907
#      pass eval_w100 indel    72    67   54 13  18     0.806   0.696   0.883  0.750  0.639  0.836 0.777
#  pipeline      eval   SNV  1159  1070 1047 23 112     0.979   0.968   0.986  0.903  0.885  0.919 0.939
#  pipeline      eval indel    48    56   41 15   7     0.732   0.604   0.830  0.854  0.728  0.928 0.788
#  pipeline eval_w100   SNV  1666  1448 1402 46 264     0.968   0.958   0.976  0.842  0.823  0.858 0.900
#  pipeline eval_w100 indel    72    73   54 19  18     0.740   0.629   0.827  0.750  0.639  0.836 0.745
####

# Step 5. Read evidence at every scored truth position
# What the caller saw: tumour and normal reads at each truth site in the recalibrated BAMs,
# counted with Mutect2's default read filters (mapping quality >= 20, base quality >= 10,
# duplicates, secondary and QC-fail reads removed, no BAQ). Counts are per DNA fragment:
# where a read pair overlaps, bcftools counts the two mates once, as Mutect2 does by
# default (--independent-mates false). Mutect2's FORMAT/DP counts both mates, so it runs
# higher (median 104 vs 78 at the called truth SNVs here).
te <- truth %>% filter(in_eval)
write.table(unique(te[, c("CHROM", "POS")]), w("truth_sites.txt"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
for (s in c("tumour", "normal")) {
    bam <- if (s == "tumour") tumour.bam else normal.bam
    sh(sprintf("bcftools mpileup -f %s -R %s -q 20 -Q 10 -B -d 100000 -a FORMAT/AD,FORMAT/DP -Ou %s 2>/dev/null | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t[%%AD\\t%%DP]\\n' > %s",
               shQuote(ref.fa), w("truth_sites.txt"), shQuote(bam), w(paste0("pileup_", s, ".tsv"))))
}
pile.alt <- function(s) {
    p <- read.tsv(w(paste0("pileup_", s, ".tsv")), c("CHROM", "POS", "PREF", "PALT", "AD", "DP"))
    p <- p[!duplicated(paste(p$CHROM, p$POS)), ]   # one line per position (indel lines follow)
    m <- match(paste(te$CHROM, te$POS), paste(p$CHROM, p$POS))
    alts <- strsplit(p$PALT[m], ","); ads <- strsplit(p$AD[m], ",")
    alt.n <- mapply(function(a, d, want) {
        if (is.null(a) || all(is.na(a))) return(0)
        i <- match(want, a); if (is.na(i)) 0 else as.numeric(d[i + 1])
    }, alts, ads, te$ALT)
    data.frame(depth = as.numeric(ifelse(is.na(m), 0, p$DP[m])), alt = alt.n)
}
pt <- pile.alt("tumour"); pn <- pile.alt("normal")
te$T_depth <- pt$depth; te$N_depth <- pn$depth
te$T_alt <- ifelse(te$type == "SNV", pt$alt, NA)   # pileup allele counts are reliable for SNVs only
te$N_alt <- ifelse(te$type == "SNV", pn$alt, NA)
m <- match(te$key, calls$key)
te$call_filter <- calls$FILTER[m]
te$called <- !is.na(m) & te$call_filter %in% "PASS"
te$status <- ifelse(te$called, "called (PASS)", ifelse(is.na(m), "not emitted", "emitted, filtered"))
cat("\nStep 5. Tumour depth at scored truth sites (median, 10th and 90th percentile)\n")
print(as.data.frame(te %>% group_by(type) %>% summarise(n = n(), median = median(T_depth),
      p10 = quantile(T_depth, 0.1), p90 = quantile(T_depth, 0.9), .groups = "drop")), row.names = FALSE)

# Check the pileup against Mutect2 at the called truth SNVs: depth and VAF should agree,
# and the observed tumour VAF should track the SEQC2 VAF from deep WGS.
chk <- te %>% filter(type == "SNV", status == "called (PASS)") %>%
    mutate(M_DP = calls$T_DP[match(key, calls$key)], M_ALT = calls$T_ALT[match(key, calls$key)])
cat(sprintf("  called truth SNVs: median depth %.0f fragments vs Mutect2 DP %.0f reads; variant fragments vs Mutect2 alt reads r = %.3f\n",
            median(chk$T_depth), median(chk$M_DP), cor(chk$T_alt, chk$M_ALT)))
vaf.ok <- te %>% filter(type == "SNV", T_depth >= 20)
cat(sprintf("  observed tumour VAF vs SEQC2 TVAF at truth SNVs with depth >= 20 (n = %d): r = %.3f\n",
            nrow(vaf.ok), cor(vaf.ok$T_alt / vaf.ok$T_depth, vaf.ok$TVAF)))
write.table(te, file.path(out.dir, "truth_sites_scored.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

####
# Step 5. Tumour depth at scored truth sites (median, 10th and 90th percentile)
#   type    n median  p10   p90
#    SNV 1159     76 29.0 185.0
#  indel   48     73 29.2 203.3
#   called truth SNVs: median depth 78 fragments vs Mutect2 DP 104 reads; variant fragments vs Mutect2 alt reads r = 0.994
#   observed tumour VAF vs SEQC2 TVAF at truth SNVs with depth >= 20 (n = 1111): r = 0.987
####

# Step 6. Why truth SNVs were missed
# Each missed SNV is either emitted and then removed by FilterMutectCalls (listed by filter),
# or never emitted, split by how many tumour fragments carry the variant.
fn <- te %>% filter(type == "SNV", status != "called (PASS)") %>%
    mutate(reason = case_when(
        status == "emitted, filtered" ~ paste("filtered:", call_filter),
        T_alt == 0 ~ "not emitted: 0 alt fragments",
        T_alt <= 2 ~ "not emitted: 1-2 alt fragments",
        TRUE ~ "not emitted: 3+ alt fragments"))
cat("\nStep 6. Missed truth SNVs in the scored region, by reason\n")
print(as.data.frame(fn %>% count(reason, sort = TRUE) %>% mutate(share = round(n / sum(n), 3))), row.names = FALSE)
cat("  missed SNVs: median tumour depth", median(fn$T_depth), "; median SEQC2 TVAF", round(median(fn$TVAF), 3),
    "(called:", round(median(te$TVAF[te$type == "SNV" & te$called]), 3), ")\n")
cont <- fn %>% filter(grepl("contamination", call_filter)) %>% mutate(POPAF = calls$POPAF[match(key, calls$key)])
cat(sprintf("  removed with the contamination filter: %d; SEQC2 TVAF 10th, 50th, 90th percentile %.3f, %.3f, %.3f; absent from gnomAD (POPAF 6) %d\n",
            nrow(cont), quantile(cont$TVAF, 0.1), median(cont$TVAF), quantile(cont$TVAF, 0.9), sum(cont$POPAF >= 6)))

# Figure 1. Every scored truth SNV by tumour depth and SEQC2 VAF, coloured by outcome.
# The dashed line is depth x VAF = 3: below it fewer than three variant fragments are
# expected. Colours are the first three categorical slots of the dataviz reference palette
# (validated all-pairs for colour-vision deficiency); counts are in the legend.
pal <- c("called (PASS)" = "#2a78d6", "emitted, filtered" = "#eb6834", "not emitted" = "#1baf7a")
ink <- "#0b0b0b"; ink2 <- "#52514e"; grid.col <- "#e6e5e1"; surface <- "#fcfcfb"
theme.bench <- theme_minimal(base_size = 11) + theme(
    plot.background = element_rect(fill = surface, colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = grid.col, linewidth = 0.3),
    axis.text = element_text(colour = ink2), axis.title = element_text(colour = ink2),
    plot.title = element_text(colour = ink, face = "bold"), plot.subtitle = element_text(colour = ink2),
    legend.position = "top", legend.justification = "left", legend.title = element_blank(),
    legend.text = element_text(colour = ink2), plot.title.position = "plot")
fig1.d <- te %>% filter(type == "SNV") %>%
    mutate(status = factor(status, levels = names(pal)), depth = pmax(T_depth, 1)) %>% arrange(status)
n.status <- table(fig1.d$status)
lab.status <- setNames(sprintf("%s (%d)", names(pal), as.integer(n.status[names(pal)])), names(pal))
# direct labels placed in open space for this data set (checked on the rendered figure)
centres <- data.frame(status = names(pal), depth = c(60, 230, 7), TVAF = c(0.37, 0.033, 0.12))
floor.note <- if (any(fig1.d$depth < 4)) "; sites below 4x drawn at 4x" else ""
fig1 <- ggplot(fig1.d, aes(depth, TVAF)) +
    geom_function(fun = function(x) 3 / x, linetype = "dashed", colour = ink2, linewidth = 0.4) +
    annotate("text", x = 420, y = 3 / 420, label = "3 variant fragments expected", colour = ink2,
             size = 3, hjust = 1, vjust = -0.6, angle = -25) +
    geom_point(aes(fill = status), shape = 21, colour = surface, stroke = 0.3, size = 2.1) +
    geom_label(data = centres, aes(label = status), fill = surface, colour = ink, linewidth = 0,
               size = 3.2, alpha = 0.85) +
    scale_fill_manual(values = pal, labels = lab.status) +
    scale_x_log10(breaks = c(5, 10, 20, 50, 100, 200, 500), limits = c(4, 700), oob = scales::squish) +
    scale_y_log10(breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1), labels = scales::percent) +
    labs(x = "Tumour depth at the site (fragments, log scale)", y = "SEQC2 tumour VAF (log scale)",
         title = "Truth SNVs in the scored region, by outcome",
         subtitle = sprintf("SEQC2 HCC1395 exome; capture targets inside the high-confidence regions (%.1f Mb)%s",
                            regions$Mb[regions$region == "eval"], floor.note)) +
    theme.bench
ggsave(file.path(out.dir, "figure1_truth_snvs_by_outcome.png"), fig1, width = 7.5, height = 5.2, dpi = 150, bg = surface)

# Recall by SEQC2 VAF, by tumour depth and by tumour variant fragments
snv <- te %>% filter(type == "SNV") %>%
    mutate(vaf_bin = cut(TVAF, c(0, 0.05, 0.1, 0.2, 0.4, 1), right = FALSE, include.lowest = TRUE,
                         labels = c("<5%", "5-10%", "10-20%", "20-40%", ">=40%")),
           depth_bin = cut(T_depth, c(-Inf, 19, 49, 99, 199, Inf),
                           labels = c("<20", "20-49", "50-99", "100-199", ">=200")),
           alt_bin = cut(T_alt, c(-Inf, 0, 2, 4, 9, Inf), labels = c("0", "1-2", "3-4", "5-9", ">=10")))
strat <- function(v) {
    d <- snv %>% group_by(level = .data[[v]]) %>% summarise(truth = n(), called = sum(called), .groups = "drop")
    ci <- t(mapply(wilson, d$called, d$truth))
    d %>% mutate(stratum = v, recall = called / truth, lo = ci[, 1], hi = ci[, 2], .before = 1)
}
recall.strata <- bind_rows(strat("vaf_bin"), strat("depth_bin"), strat("alt_bin")) %>%
    mutate(level = as.character(level))
cat("\n  Recall of truth SNVs by stratum\n")
print(as.data.frame(recall.strata %>% mutate(across(c(recall, lo, hi), ~ round(.x, 3)))), row.names = FALSE)
write.table(recall.strata, file.path(out.dir, "recall_by_stratum.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Figure 2. Recall by SEQC2 VAF, with 95% Wilson intervals
fig2.d <- recall.strata %>% filter(stratum == "vaf_bin") %>%
    mutate(level = factor(level, levels = levels(snv$vaf_bin)))
fig2 <- ggplot(fig2.d, aes(level, recall)) +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = pal[["called (PASS)"]], linewidth = 0.7) +
    geom_point(shape = 21, fill = pal[["called (PASS)"]], colour = surface, stroke = 0.6, size = 3.2) +
    geom_text(aes(y = hi, label = sprintf("%d of %d", called, truth)), vjust = -0.8, colour = ink2, size = 3.2) +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2), labels = scales::percent) +
    labs(x = "SEQC2 tumour VAF (pooled deep WGS)", y = "Recall, PASS calls",
         title = "Recall of truth SNVs by variant allele fraction",
         subtitle = "Scored region; bars are 95% Wilson intervals, labels are called of truth") +
    theme.bench + theme(panel.grid.major.x = element_blank())
ggsave(file.path(out.dir, "figure2_recall_by_vaf.png"), fig2, width = 6.5, height = 4.2, dpi = 150, bg = surface)

####
# Step 6. Missed truth SNVs in the scored region, by reason
#                                             reason  n share
#                            filtered: contamination 39 0.348
#                     not emitted: 1-2 alt fragments 22 0.196
#                       not emitted: 0 alt fragments  9 0.080
#              filtered: contamination;weak_evidence  8 0.071
#                              filtered: strand_bias  8 0.071
#                filtered: contamination;strand_bias  6 0.054
#                      not emitted: 3+ alt fragments  6 0.054
#                              filtered: orientation  5 0.045
#                         filtered: clustered_events  2 0.018
#            filtered: clustered_events;multiallelic  1 0.009
#                filtered: contamination;orientation  1 0.009
#    filtered: contamination;orientation;strand_bias  1 0.009
#  filtered: contamination;orientation;weak_evidence  1 0.009
#  filtered: contamination;strand_bias;weak_evidence  1 0.009
#                filtered: strand_bias;weak_evidence  1 0.009
#                            filtered: weak_evidence  1 0.009
#   missed SNVs: median tumour depth 51 ; median SEQC2 TVAF 0.069 (called: 0.373 )
#   removed with the contamination filter: 57; SEQC2 TVAF 10th, 50th, 90th percentile 0.038, 0.058, 0.104; absent from gnomAD (POPAF 6) 51
#
#   Recall of truth SNVs by stratum
#    stratum recall    lo    hi   level truth called
#    vaf_bin  0.237 0.130 0.392     <5%    38      9
#    vaf_bin  0.607 0.518 0.689   5-10%   122     74
#    vaf_bin  0.886 0.829 0.926  10-20%   167    148
#    vaf_bin  0.972 0.947 0.985  20-40%   320    311
#    vaf_bin  0.986 0.972 0.993   >=40%   512    505
#  depth_bin  0.708 0.568 0.818     <20    48     34
#  depth_bin  0.863 0.818 0.899   20-49   278    240
#  depth_bin  0.924 0.895 0.945   50-99   432    399
#  depth_bin  0.938 0.906 0.960 100-199   308    289
#  depth_bin  0.914 0.839 0.956   >=200    93     85
#    alt_bin  0.000 0.000 0.299       0     9      0
#    alt_bin  0.171 0.085 0.313     1-2    41      7
#    alt_bin  0.400 0.276 0.538     3-4    50     20
#    alt_bin  0.807 0.735 0.863     5-9   145    117
#    alt_bin  0.988 0.979 0.993    >=10   914    903
####

# Step 7. False positives
# PASS calls in the scored region that are not in the reference set: substitution, read
# support, read-orientation split of the alt reads, whether they came from a split
# multi-base call, and whether any SEQC2 caller reported them (the superset carries every
# call any SEQC2 pipeline made, graded HighConf / MedConf / LowConf / Unclassified).
fp <- calls %>% filter(keep_pass, in_eval, !key %in% truth$key[truth$in_eval])
fp.pipe <- calls %>% filter(keep_pipeline, !keep_pass, in_eval, !key %in% truth$key[truth$in_eval])
look <- bind_rows(fp, fp.pipe)
write.table(unique(look[, c("CHROM", "POS")]), w("fp_sites.txt"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
for (t in c("snv", "indel")) {
    ss <- if (t == "snv") superset.snv else superset.indel
    sh(sprintf("bcftools view -R %s %s -Ou | %s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%FILTER\\n' > %s",
               w("fp_sites.txt"), shQuote(ss), norm, w(paste0("superset_fp_", t, ".tsv"))))
}
sup <- bind_rows(read.tsv(w("superset_fp_snv.tsv"), c("CHROM", "POS", "REF", "ALT", "SUPERSET")),
                 read.tsv(w("superset_fp_indel.tsv"), c("CHROM", "POS", "REF", "ALT", "SUPERSET")))
sup$key <- key(sup)
comp <- c(A = "T", C = "G", G = "C", T = "A")
fp.table <- function(d) {
    d %>% mutate(
        superset = ifelse(key %in% sup$key, sub("^PASS;", "", sup$SUPERSET[match(key, sup$key)]), "absent"),
        from_mnp = !is.na(ORIGINAL) & sapply(strsplit(ORIGINAL, "|", fixed = TRUE), function(o)
            length(o) >= 4 && nchar(o[3]) > 1 && nchar(o[3]) == nchar(o[4])),
        # substitution on the pyrimidine strand (C>A also stands for G>T)
        sub = ifelse(type == "SNV", ifelse(REF %in% c("C", "T"), paste0(REF, ">", ALT),
                                           paste0(comp[REF], ">", comp[ALT])), "indel"),
        alt_F1R2 = alt.of(T_F1R2), alt_F2R1 = alt.of(T_F2R1),
        orient_frac = pmax(alt_F1R2, alt_F2R1) / (alt_F1R2 + alt_F2R1)) %>%
        select(CHROM, POS, REF, ALT, FILTER, type, sub, from_mnp, T_AF, T_DP, T_ALT, N_ALT, TLOD, POPAF,
               alt_F1R2, alt_F2R1, orient_frac, superset)
}
fp.out <- fp.table(fp)
cat("\nStep 7. PASS false positives in the scored region:", nrow(fp.out), "\n")
short <- function(x) ifelse(nchar(x) > 8, paste0(substr(x, 1, 6), "..", nchar(x)), x)   # long alleles: first 6 bases and length
print(as.data.frame(fp.out %>% mutate(REF = short(REF), ALT = short(ALT)) %>%
      select(CHROM, POS, REF, ALT, sub, from_mnp, T_AF, T_DP, T_ALT, N_ALT, TLOD, alt_F1R2, alt_F2R1, superset)),
      row.names = FALSE, digits = 3)
tp.sub <- calls %>% filter(keep_pass, in_eval, type == "SNV", key %in% truth$key) %>%
    mutate(sub = ifelse(REF %in% c("C", "T"), paste0(REF, ">", ALT), paste0(comp[REF], ">", comp[ALT])))
cat(sprintf("  C>A (G>T) share: %d of %d false-positive SNVs vs %.3f of true-positive SNVs\n",
            sum(fp.out$sub == "C>A"), sum(fp.out$type == "SNV"), mean(tp.sub$sub == "C>A")))
cat("  extra false positives admitted by the pipeline keep rule (normal_artifact):", nrow(fp.pipe), "\n")
print(as.data.frame(fp.table(fp.pipe) %>% mutate(REF = short(REF), ALT = short(ALT)) %>%
      select(CHROM, POS, REF, ALT, type, T_AF, T_ALT, N_ALT, superset)),
      row.names = FALSE, digits = 3)
write.table(bind_rows(fp.out %>% mutate(callset = "pass"), fp.table(fp.pipe) %>% mutate(callset = "pipeline only")),
            file.path(out.dir, "false_positives.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

####
# Step 7. PASS false positives in the scored region: 20
#  CHROM       POS        REF        ALT   sub from_mnp  T_AF T_DP T_ALT N_ALT   TLOD alt_F1R2 alt_F2R1 superset
#   chr1 167425530          C CTGTCT..45 indel    FALSE 0.259   41     7     0  23.24        3        3   absent
#   chr2 109129817          C          T   C>T    FALSE 0.148   31     5     0   8.98        1        2   absent
#   chr2 197500318      CCGCT          C indel    FALSE 0.427   86    35     0 131.59       14       12   absent
#   chr2 197500324 GGGCGA..12          G indel    FALSE 0.391   95    35     0 129.93       14       12   absent
#   chr5 127417888          C          A   C>A    FALSE 0.090   36     4     0   5.66        1        1   absent
#   chr7 129611672          G          T   C>A    FALSE 0.248   13     3     0   6.39        1        1   absent
#   chr8  64581002 AGGGGG..44          A indel    FALSE 0.266   80    19     0  68.46       10        5   absent
#   chr9  41644642         CG          C indel    FALSE 0.273   13     3     0   7.41        1        1   absent
#   chr9 136197714          C          T   C>T    FALSE 0.060   78     5     0   7.08        1        2   absent
#  chr15  43220669          T          C   T>C    FALSE 0.250   10     2     0   6.31        0        1   absent
#  chr15  74409355          G          A   C>T    FALSE 0.088   35     3     0   8.00        1        1   absent
#  chr16  19551340        TTA          T indel    FALSE 0.251  526   129     0 476.96       27       59   absent
#  chr16  19551343 CACAGA..24          C indel    FALSE 0.219  635   129     0 470.26       27       59   absent
#  chr16  28470455          G          T   C>A    FALSE 0.182   13     2     0   7.05        1        0   absent
#  chr19   6833167         TA          T indel    FALSE 0.305   12     4     0   9.64        0        3   absent
#  chr20  17693692        CTG          C indel    FALSE 0.293   48    14     0  51.17        2        5   absent
#  chr20  17693696 GAGGAA..42          G indel    FALSE 0.267   57    14     0  50.60        2        4   absent
#  chr20  30388434      TATCA          T indel    FALSE 0.146   60     7     0  20.22        4        3   absent
#  chr20  59840349          C          T   C>T     TRUE 0.171   34     5     0  12.15        1        3   absent
#  chr20  59840350          C          T   C>T     TRUE 0.171   34     5     0  12.15        1        3   absent
#   C>A (G>T) share: 3 of 9 false-positive SNVs vs 0.130 of true-positive SNVs
#   extra false positives admitted by the pipeline keep rule (normal_artifact): 18
#  CHROM       POS         REF ALT  type  T_AF T_ALT N_ALT superset
#   chr1  16627332           G  GC indel 0.208    54    31   absent
#   chr1  17393975  GATGAG..93   G indel 0.066    38    16   absent
#   chr3    406414          TC   T indel 0.139    23    16   absent
#   chr3 197621775           T   C   SNV 0.135    14    11   absent
#   chr4 189983293           C   T   SNV 0.137    23    13   absent
#   chr5 173151565           G   T   SNV 0.234     8     2   absent
#   chr8 100709464           G   A   SNV 0.122    21     8   absent
#   chr9  62801804           G   A   SNV 0.152    28    23   absent
#  chr11  68021828           T   C   SNV 0.042    17    17   absent
#  chr12  31112889           A   G   SNV 0.202    18     9   absent
#  chr12  40486899           C   T   SNV 0.078    51    80   absent
#  chr13  25096984           C   T   SNV 0.310    53    21   absent
#  chr15  39589305           G   A   SNV 0.081     4     1   absent
#  chr17  21416665           A   G   SNV 0.229    27    15   absent
#  chr21  10473366           A   T   SNV 0.084    42    60   absent
#  chr21  10473402           T   C   SNV 0.076    50    59   absent
#  chr21  10491351           T   C   SNV 0.062    13     9   absent
#  chr21  44591681 GTGCTG..127   G indel 0.025    20    30   absent
####

# Step 8. Indel representation check
# An unmatched indel call within 10 bp of an unmatched truth indel may be the same event
# written differently; counted here, not credited in Step 4.
fn.ind <- te %>% filter(type == "indel", status != "called (PASS)")
fp.ind <- fp.out %>% filter(type == "indel")
near <- sapply(seq_len(nrow(fp.ind)), function(i)
    any(fn.ind$CHROM == fp.ind$CHROM[i] & abs(as.numeric(fn.ind$POS) - as.numeric(fp.ind$POS[i])) <= 10))
cat("\nStep 8. Indels: missed", nrow(fn.ind), "; false positives", nrow(fp.ind),
    "; false positives within 10 bp of a missed truth indel:", sum(near), "\n")
print(as.data.frame(fn.ind %>% mutate(REF = short(REF), ALT = short(ALT)) %>%
      select(CHROM, POS, REF, ALT, TVAF, T_depth, status, call_filter)), row.names = FALSE, digits = 3)

####
# Step 8. Indels: missed 7 ; false positives 11 ; false positives within 10 bp of a missed truth indel: 0
#  CHROM       POS        REF ALT  TVAF T_depth            status                 call_filter
#   chr7  44020406 GAAGAC..16   G 0.383       8       not emitted                        <NA>
#   chr7 108106624 GTCAAA..23   G 0.082      33       not emitted                        <NA>
#   chr9   6554804 AAGCAA..32   A 0.065      24 emitted, filtered contamination;weak_evidence
#  chr12 122508614 ATTTAT..27   A 0.034      31 emitted, filtered               contamination
#  chr15  52112320         CG   C 0.067      48       not emitted                        <NA>
#  chr18  36355492 CTCCAT..25   C 0.029      90 emitted, filtered contamination;weak_evidence
#  chr19  15114940 CCATCT..28   C 0.044      49 emitted, filtered               contamination
####

# Step 9. PASS calls the truth set cannot score
# Calls outside the high-confidence regions are neither right nor wrong by this benchmark.
pass.all <- calls %>% filter(keep_pass)
cat("\nStep 9. Where the", nrow(pass.all), "PASS calls fall\n")
print(as.data.frame(pass.all %>% mutate(where = case_when(
    in_eval ~ "scored (targets n HC)",
    in_targets ~ "targets, outside HC",
    in_eval_w100 ~ "padding, inside HC",
    in_targets_w100 ~ "padding, outside HC",
    TRUE ~ "outside the padded targets")) %>% count(where, type) %>% arrange(where, type)), row.names = FALSE)
# Density per Mb, and whether any SEQC2 caller reported the PASS calls in targets outside HC
# (the superset is genome-wide, so it covers these regions too).
dens <- function(r) sum(pass.all[[paste0("in_", r)]]) / regions$Mb[regions$region == r]
out.hc <- pass.all %>% filter(in_targets, !in_eval)
write.table(unique(out.hc[, c("CHROM", "POS")]), w("outhc_sites.txt"), sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
for (t in c("snv", "indel")) {
    ss <- if (t == "snv") superset.snv else superset.indel
    sh(sprintf("bcftools view -R %s %s -Ou | %s | bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%FILTER\\n' > %s",
               w("outhc_sites.txt"), shQuote(ss), norm, w(paste0("superset_outhc_", t, ".tsv"))))
}
sup2 <- bind_rows(read.tsv(w("superset_outhc_snv.tsv"), c("CHROM", "POS", "REF", "ALT", "SUPERSET")),
                  read.tsv(w("superset_outhc_indel.tsv"), c("CHROM", "POS", "REF", "ALT", "SUPERSET")))
out.hc$superset <- ifelse(key(out.hc) %in% key(sup2), sub("^PASS;", "", sup2$SUPERSET[match(key(out.hc), key(sup2))]), "absent")
targets.out.mb <- regions$Mb[regions$region == "targets"] - regions$Mb[regions$region == "eval"]
cat(sprintf("  PASS calls per Mb: scored region %.1f; targets outside HC %.1f (%.2f Mb)\n",
            dens("eval"), nrow(out.hc) / targets.out.mb, targets.out.mb))
cat("  PASS calls in targets outside HC, by SEQC2 superset grade (HighConf in any tier vs the rest)\n")
print(as.data.frame(out.hc %>% mutate(grade = case_when(grepl("HighConf", superset) ~ "HighConf",
                                                        superset == "absent" ~ "absent",
                                                        TRUE ~ "MedConf, LowConf or Unclassified")) %>%
      count(type, grade) %>% group_by(type) %>% mutate(share = round(n / sum(n), 3)) %>% ungroup()), row.names = FALSE)
# Most of them sit where the normal line carries a single parental copy: common germline
# variants (gnomAD AF >= 1%, POPAF < 2) with no support in the normal, so they look somatic.
cat("  The same calls by chromosome: common gnomAD variants and normal support\n")
loh <- bind_rows(out.hc %>% mutate(group = ifelse(CHROM %in% c("chr6", "chr16", "chrX"), CHROM, "other")),
                 pass.all %>% filter(in_eval) %>% mutate(group = "scored region")) %>%
    group_by(group) %>%
    summarise(calls = n(), gnomad_common = sum(POPAF < 2, na.rm = TRUE),
              normal_alt_zero = sum(N_ALT == 0 & POPAF < 2, na.rm = TRUE),
              mean_tumour_vaf = round(mean(T_AF), 3),
              span_Mb = sprintf("%.1f-%.1f", min(as.numeric(POS[POPAF < 2])) / 1e6, max(as.numeric(POS[POPAF < 2])) / 1e6),
              .groups = "drop")
print(as.data.frame(loh), row.names = FALSE)

####
# Step 9. Where the 3648 PASS calls fall
#                  where  type    n
#     padding, inside HC   SNV  371
#     padding, inside HC indel   15
#    padding, outside HC   SNV  451
#    padding, outside HC indel   42
#  scored (targets n HC)   SNV 1056
#  scored (targets n HC) indel   52
#    targets, outside HC   SNV 1516
#    targets, outside HC indel  145
#   PASS calls per Mb: scored region 13.7; targets outside HC 169.0 (9.83 Mb)
#   PASS calls in targets outside HC, by SEQC2 superset grade (HighConf in any tier vs the rest)
#   type                            grade    n share
#    SNV                         HighConf 1435 0.947
#    SNV MedConf, LowConf or Unclassified   68 0.045
#    SNV                           absent   13 0.009
#  indel                         HighConf   73 0.503
#  indel MedConf, LowConf or Unclassified   54 0.372
#  indel                           absent   18 0.124
#   The same calls by chromosome: common gnomAD variants and normal support
#          group calls gnomad_common normal_alt_zero mean_tumour_vaf   span_Mb
#          chr16   446           390             378           0.522 46.6-89.9
#           chr6   817           724             701           0.680  0.4-56.6
#           chrX   339           274             260           0.921 0.3-155.9
#          other    59             2               0           0.409 43.7-54.2
#  scored region  1108             3               3           0.449 22.9-81.9
####

# Step 10. Independent recount
# bcftools isec on the same normalised records must reproduce the Step 4 counts exactly.
recount <- bind_rows(lapply(c("SNV", "indel"), function(t) {
    tv <- if (t == "SNV") truth.snv else truth.indel
    vt <- if (t == "SNV") "snps" else "indels"
    sh(sprintf("bcftools view -T %s %s -Ou | %s | bcftools view -v %s -Oz -o %s && bcftools index -f %s",
               w("eval.bed"), shQuote(tv), norm, vt, w("isec_truth.vcf.gz"), w("isec_truth.vcf.gz")))
    sh(sprintf("bcftools view -f PASS %s -Ou | %s | bcftools view -T %s -v %s -Oz -o %s && bcftools index -f %s",
               shQuote(calls.vcf), norm, w("eval.bed"), vt, w("isec_calls.vcf.gz"), w("isec_calls.vcf.gz")))
    sh(sprintf("bcftools isec -c none -n=2 -w1 %s %s -Ou | bcftools view -H | wc -l > %s",
               w("isec_calls.vcf.gz"), w("isec_truth.vcf.gz"), w("isec_tp.txt")))
    tp <- as.numeric(readLines(w("isec_tp.txt")))
    nc <- as.numeric(system(sprintf("bcftools view -H %s | wc -l", w("isec_calls.vcf.gz")), intern = TRUE))
    nt <- as.numeric(system(sprintf("bcftools view -H %s | wc -l", w("isec_truth.vcf.gz")), intern = TRUE))
    s4 <- scores %>% filter(callset == "pass", region == "eval", type == t)
    data.frame(type = t, TP = tp, FP = nc - tp, FN = nt - tp, matches_step4 = tp == s4$TP & nc - tp == s4$FP & nt - tp == s4$FN)
}))
cat("\nStep 10. bcftools isec recount, scored region, PASS\n")
print(recount, row.names = FALSE)
if (!all(recount$matches_step4)) stop("bcftools isec disagrees with the Step 4 counts")

####
# Step 10. bcftools isec recount, scored region, PASS
#   type   TP FP  FN matches_step4
#    SNV 1047  9 112          TRUE
#  indel   41 11   7          TRUE
####
