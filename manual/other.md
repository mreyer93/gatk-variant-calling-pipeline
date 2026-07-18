## Other things found in this repository
This repo also carries historical, project-specific analysis code that isn't part of the
generalized pipeline described in the rest of this manual - things like one-off investigations
into a specific cohort, RNA-seq analyses for specific experiments, and superseded execution
models (e.g. the old AWS/S3-based setup). All of that now lives under [`archive/`](../archive/README.md),
which explains what's in each subdirectory. Nothing under `archive/` is maintained, documented, or
guaranteed to run - it's kept only as reference.

`vcf_clustering_standalone` (in `archive/exploratory_analysis/`) has ad hoc scripts to cluster
samples by variant calls. Ideally this would be done with germline data, generalized, and run
automatically on each new batch of data as a proper Snakemake step - it's a reasonable follow-on
enhancement for this pipeline if you need it again.
