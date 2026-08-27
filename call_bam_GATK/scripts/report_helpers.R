# Helper functions for the somatic variant-calling report.
#
# Design notes:
#  - Every loader returns NULL (not an error) when its inputs are absent, so the report
#    degrades gracefully. The pipeline can legitimately be run in several modes:
#    full annotation, skip_annotation (see config_call_bam_GATK_local.yaml), no matched
#    normals, or CADD skipped - the report adapts rather than failing.
#  - Nothing here assumes a particular sample-naming convention; sample names come from
#    the metadata file and output filenames.

suppressPackageStartupMessages({
    library(data.table)
})

# Read a VCF-like table (the pipeline's annotated outputs are TSVs with a header row;
# raw GATK VCFs have ## meta-lines then a #CHROM header). Returns NULL on missing/empty.
read_variant_table <- function(path) {
    if (!file.exists(path)) return(NULL)
    if (file.info(path)$size == 0) return(NULL)
    out <- tryCatch({
        first <- readLines(path, n = 1, warn = FALSE)
        if (length(first) == 0) return(NULL)
        if (startsWith(first, "##")) {
            # a real VCF - let fread find the #CHROM header line
            dt <- data.table::fread(path, skip = "#CHROM", sep = "\t",
                                    quote = "", showProgress = FALSE)
        } else {
            dt <- data.table::fread(path, sep = "\t", quote = "", showProgress = FALSE)
        }
        if (nrow(dt) == 0) return(NULL)
        setnames(dt, 1, sub("^#", "", names(dt)[1]))
        dt
    }, error = function(e) {
        warning("Could not read ", path, ": ", conditionMessage(e))
        NULL
    })
    out
}

# Load per-sample annotated variant tables produced by combine_annotation.R.
# Returns one data.table with a `sample` column, or NULL if the directory has nothing.
load_annotated_variants <- function(annot_dir, suffix = "_filtered_annotated.vcf") {
    if (!dir.exists(annot_dir)) return(NULL)
    files <- list.files(annot_dir, pattern = paste0(suffix, "$"), full.names = TRUE)
    if (length(files) == 0) return(NULL)
    parts <- lapply(files, function(f) {
        dt <- read_variant_table(f)
        if (is.null(dt)) return(NULL)
        dt[, sample := sub(suffix, "", basename(f))]
        # AF_/DP_ columns are named per sample; normalise them so samples can be stacked
        af_col <- grep("^AF_", names(dt), value = TRUE)
        dp_col <- grep("^DP_", names(dt), value = TRUE)
        if (length(af_col) == 1) dt[, AF := suppressWarnings(as.numeric(sub(",.*$", "", get(af_col))))]
        if (length(dp_col) == 1) dt[, DP := suppressWarnings(as.numeric(sub(",.*$", "", get(dp_col))))]
        keep <- intersect(c("sample", "CHROM", "POS", "REF", "ALT", "FILTER", "hugoSymbol",
                            "variantClassification", "variantType", "proteinChange",
                            "CADD_phred", "AF", "DP"), names(dt))
        dt[, ..keep]
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0) return(NULL)
    data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
}

# Load the per-patient ranked "new score" tables (combine_samples.R / combine_samples_redo.R)
load_new_scores <- function(scores_dir) {
    if (!dir.exists(scores_dir)) return(NULL)
    files <- list.files(scores_dir, pattern = "\\.vcf$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    parts <- lapply(files, function(f) {
        dt <- read_variant_table(f)
        if (is.null(dt)) return(NULL)
        dt[, patient := sub("\\.vcf$", "", basename(f))]
        dt
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0) return(NULL)
    data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
}

# Count variants straight from called VCFs. Used when annotation was skipped, so the
# report still has something quantitative to say about the calls.
load_vcf_counts <- function(vcf_dir, suffix = "_filtered.vcf") {
    if (!dir.exists(vcf_dir)) return(NULL)
    files <- list.files(vcf_dir, pattern = paste0(suffix, "$"), full.names = TRUE)
    if (length(files) == 0) return(NULL)
    rows <- lapply(files, function(f) {
        dt <- read_variant_table(f)
        sample <- sub(suffix, "", basename(f))
        if (is.null(dt)) {
            return(data.table(sample = sample, n_variants = 0L, n_pass = 0L))
        }
        n_pass <- if ("FILTER" %in% names(dt)) sum(dt$FILTER == "PASS", na.rm = TRUE) else NA_integer_
        data.table(sample = sample, n_variants = nrow(dt), n_pass = as.integer(n_pass))
    })
    data.table::rbindlist(rows)
}

# GATK CalculateContamination output: one table per tumour sample, columns
# level / contamination / error. Returns NULL if the step didn't run.
load_contamination <- function(contam_dir) {
    if (!dir.exists(contam_dir)) return(NULL)
    files <- list.files(contam_dir, pattern = "_contam\\.table$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    rows <- lapply(files, function(f) {
        dt <- tryCatch(data.table::fread(f, showProgress = FALSE), error = function(e) NULL)
        if (is.null(dt) || nrow(dt) == 0) return(NULL)
        data.table(sample = sub("_contam\\.table$", "", basename(f)),
                   contamination = suppressWarnings(as.numeric(dt[[2]][1])),
                   error = suppressWarnings(as.numeric(dt[[3]][1])))
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)
    data.table::rbindlist(rows)
}

# Joined smoothed-depth table (one column per sample after the annotation columns)
load_depth_summary <- function(depth_file) {
    if (!file.exists(depth_file)) return(NULL)
    dt <- tryCatch(data.table::fread(depth_file, showProgress = FALSE),
                   error = function(e) NULL)
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    meta_cols <- intersect(c("name", "chr", "start", "end", "strand", "exon",
                             "length", "gene"), names(dt))
    sample_cols <- setdiff(names(dt), meta_cols)
    if (length(sample_cols) == 0) return(NULL)
    rows <- lapply(sample_cols, function(s) {
        v <- suppressWarnings(as.numeric(dt[[s]]))
        v <- v[is.finite(v)]
        if (length(v) == 0) return(NULL)
        data.table(sample = s,
                   mean_depth = mean(v),
                   median_depth = stats::median(v),
                   frac_below_50x = mean(v < 50))
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)
    data.table::rbindlist(rows)
}

# FACETS purity/ploidy estimates, one small text file per patient-timepoint
load_facets_purity <- function(facets_dir) {
    if (!dir.exists(facets_dir)) return(NULL)
    files <- list.files(facets_dir, pattern = "_purity\\.txt$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    rows <- lapply(files, function(f) {
        txt <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
        if (length(txt) == 0) return(NULL)
        nums <- suppressWarnings(as.numeric(regmatches(txt, regexpr("[0-9.]+", txt))))
        data.table(patient_tp = sub("_purity\\.txt$", "", basename(f)),
                   value = paste(txt, collapse = "; "),
                   first_numeric = if (length(nums)) nums[1] else NA_real_)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)
    data.table::rbindlist(rows)
}

# Consistent, colourblind-safe-ish palette and a shared theme so every figure in the
# report looks like it came from the same document.
report_theme <- function() {
    ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
            plot.title = ggplot2::element_text(face = "bold", size = 12)
        )
}

REPORT_PALETTE <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
                    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD")

# Small helper for the report: emit a short italic note when a section has no data,
# instead of erroring or silently rendering an empty plot.
note_missing <- function(what) {
    cat("\n*", what, "was not available for this run - this section is omitted.*\n\n")
}
