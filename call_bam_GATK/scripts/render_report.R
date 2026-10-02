# Renders an R Markdown report from a Snakemake rule.
#
# Kept separate from the .Rmd itself so that (a) HTML and PDF can be rendered by
# independent rules - a broken LaTeX toolchain then costs you the PDF but not the HTML -
# and (b) each format renders entirely inside its own scratch directory (figures, LaTeX
# files and the report itself), and only the finished report is moved into place.
# Snakemake runs the two rules concurrently; rendered straight into the shared output
# directory, both wrote figures into the same <name>_files/ folder, which rmarkdown
# deletes when a render finishes.

log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
    sink(type = "message")
    sink(type = "output")
    close(log_con)
}, add = TRUE)

rmd_in   <- normalizePath(snakemake@params[["rmd"]], mustWork = TRUE)
out_file <- snakemake@output[[1]]
fmt      <- snakemake@params[["format"]]

out_dir <- normalizePath(dirname(out_file), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

# per-format scratch space, removed once the report is in place (kept after a failure, for
# debugging). Not via on.exit(): Snakemake runs this file at top level, where on.exit()
# never fires.
inter_dir <- file.path(out_dir, paste0(".render_", fmt))
unlink(inter_dir, recursive = TRUE)
dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)

report_params <- list(
    project_name    = snakemake@params[["project_name"]],
    outdir          = normalizePath(snakemake@params[["outdir"]], mustWork = FALSE),
    metadata        = snakemake@params[["metadata"]],
    scripts_dir     = normalizePath(snakemake@params[["scripts_dir"]], mustWork = TRUE),
    genome_version  = snakemake@params[["genome_version"]],
    filter_min_af   = snakemake@params[["filter_min_af"]],
    filter_min_dp   = snakemake@params[["filter_min_dp"]],
    filter_min_cadd = snakemake@params[["filter_min_cadd"]],
    filter_gatk     = paste(snakemake@params[["filter_gatk"]], collapse = ", "),
    top_n           = snakemake@params[["top_n"]]
)

message("Rendering ", rmd_in, " as ", fmt, " -> ", out_file)

rmarkdown::render(
    input             = rmd_in,
    output_format     = fmt,
    output_file       = basename(out_file),
    output_dir        = inter_dir,
    intermediates_dir = inter_dir,
    knit_root_dir     = getwd(),
    params            = report_params,
    envir             = new.env(parent = globalenv()),
    quiet             = FALSE
)

rendered <- file.path(inter_dir, basename(out_file))
if (!file.exists(rendered)) {
    stop("rmarkdown::render finished but ", rendered, " was not created")
}
if (!file.rename(rendered, out_file)) {
    stop("could not move ", rendered, " to ", out_file)
}
unlink(inter_dir, recursive = TRUE)
message("Wrote ", out_file)
