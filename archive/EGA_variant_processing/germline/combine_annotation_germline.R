# combine annotations with a filtered VCF file to produce a final report

# if launched from snakemake set args here
if (exists('snakemake')) {
    vcf.input.f <- snakemake@input[['vcf_funcotator']]
    cadd.input.f <- snakemake@input[['vcf_cadd']]
    vcf.output.unfilt.f <- snakemake@output[['vcf_unfilt']]
    vcf.output.filt.f <- snakemake@output[['vcf_filt']]
    vaf.min <- snakemake@params[['vaf_min']]
    depth.min <- snakemake@params[['depth_min']]
} else{
    # parse args from command line
    args = commandArgs(trailingOnly=TRUE)
    if (length(args)!=6) {
        stop("Six arguments must be supplied: funcotator_input_decomposed decomposed_GT_input cadd_input metadata_input unfiltered_output filtered_output", call.=FALSE)
    }
    vcf.input.f <- args[1]
    gt.input.f <- args[2]
    cadd.input.f <- args[3]
    metadata.f <- args[4]
    vcf.output.unfilt.f <- args[5]
    vcf.output.filt.f <- args[6]
}

if(F){
# test with EGA chrX
vcf.input.f <- '~/local_data/EGA_data/06_GDB/chr22_decomposed.vcf'
cadd.input.f <- '~/local_data/EGA_data/06_GDB/chr22_CADD_simple.tsv'
vcf.output.unfilt.f <- '~/local_data/EGA_data/06_GDB/chr22_combined_unfiltered.vcf'
vcf.output.filt.f <- '~/local_data/EGA_data/06_GDB/chr22_combined_filtered.vcf'
gt.input.f <- '~/local_data/EGA_data/06_GDB/chr22_decomposed.GT.FORMAT'
metadata.f <- '~/local_data/EGA_data/bam_metadata_uniq_samp.tsv'
}

# read metadata
metadata <- read.table(metadata.f, sep='\t', quote='', header=T, comment.char = '')
# remove samples withe "exclude.pair.filter" in the metadata as they are poor quality samples
# that don't cluster with the rest
remove.samples <- metadata$sample[metadata$exclude.pair.filter]
# read once to get position of #CHROM
vcf <- read.table(vcf.input.f, sep='\t', quote='', comment.char = '', header = F, fill=T, nrows = 10000)
start.line <- which(vcf$V1=='#CHROM')
vcf <- read.table(vcf.input.f, sep='\t', skip = (start.line-1), quote='', comment.char = '', header = T, fill=T, check.names = F)
colnames(vcf)[1] <- 'CHROM'
gt <- read.table(gt.input.f, sep='\t', quote='', header=T, check.names = F)
gt.info  <- gt[,1:2]
gt  <- gt[,3:ncol(gt)]

sort(unique(unlist(apply(gt,2, unique))))
# we have all sorts of genotypes here
# "./." "./1" ".|." ".|1" "0/." "0/0" "0/1" "0|." "0|0" "0|1" "1/." "1/1" "1|." "1|1"
# I guess I set any call with a dot to NA, as we don't know if its haplo or homo
na.sets <- c("./.","./1",".|.",".|1","0/.","0|.","1/.","1|.")
gtmat <- matrix(NA, nrow=nrow(gt), ncol=ncol(gt), dimnames = dimnames(gt))
gtmat[gt=='0/0'] <- 0
gtmat[gt=='0|0'] <- 0
gtmat[gt=='0/1'] <- 1
gtmat[gt=='0|1'] <- 1
gtmat[gt=='1/1'] <- 2
gtmat[gt=='1|1'] <- 2

missing.variants <- apply(gtmat, 1, function(x) sum(is.na(x)))
# hist(missing.variants/ ncol(gtmat))
missing.samples <- apply(gtmat, 2, function(x) sum(is.na(x)))
# hist(missing.samples)

# remove any samples with > 50% missing data
remove.samples <- c(remove.samples, names(missing.samples[missing.samples > (nrow(gtmat) * 0.5)]))
# remove problem variants with at least 80% missing samples
remove.variants <- names(missing.variants[missing.variants > (ncol(gtmat) * 0.8)])

# filter both vcf and gtmat
gtmat <- gtmat[!(rownames(gtmat) %in% remove.variants), ]
gtmat <- gtmat[,!(colnames(gtmat) %in% remove.samples)]
vcf <- vcf[!(rownames(vcf) %in% remove.variants), ]
vcf <- vcf[,!(colnames(vcf) %in% remove.samples)]
metadata <- metadata[!metadata$sample %in% remove.samples,]
ncol(gtmat) == ncol(vcf)-9
ncol(gtmat) == nrow(metadata)
nrow(gtmat) == nrow(vcf)

# now what are the numbers
missing.variants <- apply(gtmat, 1, function(x) sum(is.na(x)))
# hist(missing.variants/ ncol(gtmat))
missing.samples <- apply(gtmat, 2, function(x) sum(is.na(x)))
# hist(missing.samples)
# variants that are perfect
full.variants <- which(missing.variants ==0 & vcf$FILTER=='PASS')

# great clustering on just a single short chromosome!
# plot(hclust(dist(t(gtmat[full.variants,]))))

# we can do some condensing down to single patients
# and eliminate NA variants this way!
gtmat.pt <- matrix(NA, nrow=nrow(gtmat), ncol=length(unique(metadata$patient)), dimnames = list(rownames(gtmat), unique(metadata$patient)))

# takes in a multi column gtmat, and aggregates
# each row according to certain rules
aggregate_patient_variants <- function(pt.gt){
    if(ncol(pt.gt)==1){
        return(as.numeric(pt.gt))
    }

    # otherwise, agg to a single col
    # if a confident genotype that's different
    a <- apply(pt.gt, 1, function(x){
        nna <- sum(is.na(x))
        if(nna == length(x)){
            return(NA)
        } else if (nna==1){
            return(max(x[!is.na(x)]))
        } else if (nna==0 & all(x==x[1])){
            return(x[1])
        } else{
            return(max(x))
        }
    })
    return(a)
}

gt.agg.list <- lapply(unique(metadata$patient), function(pt){
    print(pt)
    pt.samples <- metadata$sample[metadata$patient==pt]
    print(pt.samples)
    pt.gt <- gtmat[,pt.samples, drop=F]
    return(aggregate_patient_variants(pt.gt))
})

gt.agg <- do.call(cbind, gt.agg.list)
colnames(gt.agg) <- unique(metadata$patient)
# plot(hclust(dist(t(gt.agg[full.variants,]))))


# make something useful out of this
# INFO string, we just want some info on the gene change
funcotation.str <- sapply(vcf$INFO, function(x) {
    a <- strsplit(x, split=';', useBytes = TRUE)[[1]]
    to.ret <- a[which(sapply(a, function(b) substr(b, 1, 11)) == 'FUNCOTATION')]
    to.ret <- gsub('FUNCOTATION=\\[','',to.ret)
    to.ret <- gsub('\\]','',to.ret)
    to.ret <- gsub('\\"','',to.ret)
    # if multiple funcotations, should just return the first
    to.ret <- strsplit(to.ret, ',')[[1]][1]
    return(to.ret)
})
names(funcotation.str) <- NULL

# also going to keep the gnomAD stuff now
# should really get this from the VCF header instead of hardcoding.....
func.colnames <- c("Gencode_34_hugoSymbol", "Gencode_34_ncbiBuild", "Gencode_34_chromosome", "Gencode_34_start", "Gencode_34_end", "Gencode_34_variantClassification", "Gencode_34_secondaryVariantClassification",
                   "Gencode_34_variantType", "Gencode_34_refAllele", "Gencode_34_tumorSeqAllele1", "Gencode_34_tumorSeqAllele2", "Gencode_34_genomeChange", "Gencode_34_annotationTranscript",
                   "Gencode_34_transcriptStrand", "Gencode_34_transcriptExon", "Gencode_34_transcriptPos", "Gencode_34_cDnaChange", "Gencode_34_codonChange", "Gencode_34_proteinChange",
                   "Gencode_34_gcContent", "Gencode_34_referenceContext", "Gencode_34_otherTranscripts", "Achilles_Top_Genes", "CGC_Name", "CGC_GeneID", "CGC_Chr", "CGC_Chr_Band", "CGC_Cancer_Somatic_Mut",
                   "CGC_Cancer_Germline_Mut", "CGC_Tumour_Types__(Somatic_Mutations)", "CGC_Tumour_Types_(Germline_Mutations)", "CGC_Cancer_Syndrome", "CGC_Tissue_Type", "CGC_Cancer_Molecular_Genetics",
                   "CGC_Mutation_Type", "CGC_Translocation_Partner", "CGC_Other_Germline_Mut", "CGC_Other_Syndrome/Disease", "ClinVar_HGMD_ID", "ClinVar_SYM", "ClinVar_TYPE", "ClinVar_ASSEMBLY", "ClinVar_rs",
                   "ClinVar_VCF_AF_ESP", "ClinVar_VCF_AF_EXAC", "ClinVar_VCF_AF_TGP", "ClinVar_VCF_ALLELEID", "ClinVar_VCF_CLNDISDB", "ClinVar_VCF_CLNDISDBINCL", "ClinVar_VCF_CLNDN", "ClinVar_VCF_CLNDNINCL",
                   "ClinVar_VCF_CLNHGVS", "ClinVar_VCF_CLNREVSTAT", "ClinVar_VCF_CLNSIG", "ClinVar_VCF_CLNSIGCONF", "ClinVar_VCF_CLNSIGINCL", "ClinVar_VCF_CLNVC", "ClinVar_VCF_CLNVCSO", "ClinVar_VCF_CLNVI",
                   "ClinVar_VCF_DBVARID", "ClinVar_VCF_GENEINFO", "ClinVar_VCF_MC", "ClinVar_VCF_ORIGIN", "ClinVar_VCF_RS", "ClinVar_VCF_SSR", "ClinVar_VCF_ID", "ClinVar_VCF_FILTER", "Cosmic_overlapping_mutations",
                   "CosmicFusion_fusion_genes", "CosmicFusion_fusion_id", "CosmicTissue_total_alterations_in_gene", "CosmicTissue_tissue_types_affected", "DNARepairGenes_Activity_linked_to_OMIM",
                   "DNARepairGenes_Chromosome_location_linked_to_NCBI_MapView", "DNARepairGenes_Accession_number_linked_to_NCBI_Entrez", "Familial_Cancer_Genes_Syndrome", "Familial_Cancer_Genes_Synonym",
                   "Familial_Cancer_Genes_Reference", "Gencode_XHGNC_hgnc_id", "Gencode_XRefSeq_mRNA_id", "Gencode_XRefSeq_prot_acc", "HGNC_HGNC_ID", "HGNC_Approved_Name", "HGNC_Status", "HGNC_Locus_Type",
                   "HGNC_Locus_Group", "HGNC_Previous_Symbols", "HGNC_Previous_Name", "HGNC_Synonyms", "HGNC_Name_Synonyms", "HGNC_Chromosome", "HGNC_Date_Modified", "HGNC_Date_Symbol_Changed",
                   "HGNC_Date_Name_Changed", "HGNC_Accession_Numbers", "HGNC_Enzyme_IDs", "HGNC_Entrez_Gene_ID", "HGNC_Ensembl_Gene_ID", "HGNC_Pubmed_IDs", "HGNC_RefSeq_IDs", "HGNC_Gene_Family_ID",
                   "HGNC_Gene_Family_Name", "HGNC_CCDS_IDs", "HGNC_Vega_ID", "HGNC_Entrez_Gene_ID(supplied_by_NCBI)", "HGNC_OMIM_ID(supplied_by_OMIM)", "HGNC_RefSeq(supplied_by_NCBI)",
                   "HGNC_UniProt_ID(supplied_by_UniProt)", "HGNC_Ensembl_ID(supplied_by_Ensembl)", "HGNC_UCSC_ID(supplied_by_UCSC)", "Oreganno_Build", "Oreganno_ID", "Oreganno_Values",
                   "Simple_Uniprot_uniprot_entry_name", "Simple_Uniprot_DrugBank", "Simple_Uniprot_alt_uniprot_accessions", "Simple_Uniprot_uniprot_accession", "Simple_Uniprot_GO_Biological_Process",
                   "Simple_Uniprot_GO_Cellular_Component", "Simple_Uniprot_GO_Molecular_Function", "dbSNP_ASP", "dbSNP_ASS", "dbSNP_CAF", "dbSNP_CDA", "dbSNP_CFL", "dbSNP_COMMON", "dbSNP_DSS", "dbSNP_G5",
                   "dbSNP_G5A", "dbSNP_GENEINFO", "dbSNP_GNO", "dbSNP_HD", "dbSNP_INT", "dbSNP_KGPhase1", "dbSNP_KGPhase3", "dbSNP_LSD", "dbSNP_MTP", "dbSNP_MUT", "dbSNP_NOC", "dbSNP_NOV", "dbSNP_NSF",
                   "dbSNP_NSM", "dbSNP_NSN", "dbSNP_OM", "dbSNP_OTH", "dbSNP_PM", "dbSNP_PMC", "dbSNP_R3", "dbSNP_R5", "dbSNP_REF", "dbSNP_RS", "dbSNP_RSPOS", "dbSNP_RV", "dbSNP_S3D", "dbSNP_SAO", "dbSNP_SLO",
                   "dbSNP_SSR", "dbSNP_SYN", "dbSNP_TOPMED", "dbSNP_TPA", "dbSNP_U3", "dbSNP_U5", "dbSNP_VC", "dbSNP_VLD", "dbSNP_VP", "dbSNP_WGT", "dbSNP_WTD", "dbSNP_dbSNPBuildID", "dbSNP_ID",
                   "dbSNP_FILTER", "gnomAD_exome_AF", "gnomAD_exome_AF_afr", "gnomAD_exome_AF_afr_female", "gnomAD_exome_AF_afr_male", "gnomAD_exome_AF_amr", "gnomAD_exome_AF_amr_female", "gnomAD_exome_AF_amr_male",
                   "gnomAD_exome_AF_asj", "gnomAD_exome_AF_asj_female", "gnomAD_exome_AF_asj_male", "gnomAD_exome_AF_eas", "gnomAD_exome_AF_eas_female", "gnomAD_exome_AF_eas_jpn", "gnomAD_exome_AF_eas_kor",
                   "gnomAD_exome_AF_eas_male", "gnomAD_exome_AF_eas_oea", "gnomAD_exome_AF_female", "gnomAD_exome_AF_fin", "gnomAD_exome_AF_fin_female", "gnomAD_exome_AF_fin_male", "gnomAD_exome_AF_male",
                   "gnomAD_exome_AF_nfe", "gnomAD_exome_AF_nfe_bgr", "gnomAD_exome_AF_nfe_est", "gnomAD_exome_AF_nfe_female", "gnomAD_exome_AF_nfe_male", "gnomAD_exome_AF_nfe_nwe", "gnomAD_exome_AF_nfe_onf",
                   "gnomAD_exome_AF_nfe_seu", "gnomAD_exome_AF_nfe_swe", "gnomAD_exome_AF_oth", "gnomAD_exome_AF_oth_female", "gnomAD_exome_AF_oth_male", "gnomAD_exome_AF_popmax", "gnomAD_exome_AF_raw",
                   "gnomAD_exome_AF_sas", "gnomAD_exome_AF_sas_female", "gnomAD_exome_AF_sas_male", "gnomAD_exome_ID", "gnomAD_exome_FILTER", "gnomAD_genome_AF", "gnomAD_genome_AF_afr", "gnomAD_genome_AF_afr_female",
                   "gnomAD_genome_AF_afr_male", "gnomAD_genome_AF_amr", "gnomAD_genome_AF_amr_female", "gnomAD_genome_AF_amr_male", "gnomAD_genome_AF_asj", "gnomAD_genome_AF_asj_female", "gnomAD_genome_AF_asj_male",
                   "gnomAD_genome_AF_eas", "gnomAD_genome_AF_eas_female", "gnomAD_genome_AF_eas_male", "gnomAD_genome_AF_female", "gnomAD_genome_AF_fin", "gnomAD_genome_AF_fin_female", "gnomAD_genome_AF_fin_male",
                   "gnomAD_genome_AF_male", "gnomAD_genome_AF_nfe", "gnomAD_genome_AF_nfe_est", "gnomAD_genome_AF_nfe_female", "gnomAD_genome_AF_nfe_male", "gnomAD_genome_AF_nfe_nwe", "gnomAD_genome_AF_nfe_onf",
                   "gnomAD_genome_AF_nfe_seu", "gnomAD_genome_AF_oth", "gnomAD_genome_AF_oth_female", "gnomAD_genome_AF_oth_male", "gnomAD_genome_AF_popmax", "gnomAD_genome_AF_raw", "gnomAD_genome_ID",
                   "gnomAD_genome_FILTER")
keep.colnames <- c('Gencode_34_hugoSymbol', 'Gencode_34_ncbiBuild', 'Gencode_34_variantClassification', 'Gencode_34_secondaryVariantClassification',
               'Gencode_34_variantType', 'Gencode_34_proteinChange', 'Gencode_34_genomeChange', 'Gencode_34_cDnaChange','Gencode_34_codonChange',
               'dbSNP_ID', 'gnomAD_exome_AF', 'gnomAD_genome_AF')

# funcotation.list <- lapply(funcotation.str, function(x) strsplit(x, split='|', fixed=T, useBytes = TRUE)[[1]][1:length(func.colnames)])
funcotation.df <- do.call(rbind, lapply(funcotation.str, function(x) strsplit(x, split='|', fixed=T, useBytes = TRUE)[[1]][1:length(func.colnames)]))
rownames(funcotation.df) <- NULL
colnames(funcotation.df) <- func.colnames[1:ncol(funcotation.df)]

# subset to the columns we want to keep
funcotation.df.keep <- funcotation.df[,keep.colnames[keep.colnames %in% colnames(funcotation.df)]]
colnames(funcotation.df.keep) <- gsub('Gencode_34_', '', colnames(funcotation.df.keep))

# add this back to the main df
vcf.comb <- data.frame(vcf[, c('CHROM', 'POS', 'REF', 'ALT', 'FILTER')], as.data.frame(funcotation.df.keep), as.data.frame(gt.agg))

# get CADD score information
cadd <- read.table(cadd.input.f, sep='\t', quote='', comment.char = '', header = T, fill=T, skip = 1)
colnames(cadd)[1] <- 'CHROM'
# unique names for cadd variants and other variants
cadd$uniq.var <- paste(gsub('chr', '', cadd$CHROM), cadd$Pos, cadd$Ref, cadd$Alt, sep='_')
# remove duplicates here because all I care about is the phred score
cadd <- cadd[!duplicated(cadd$uniq.var), ]
rownames(cadd) <- cadd$uniq.var
vcf.uniq.var <- paste(gsub('chr', '', vcf.comb$CHROM), vcf.comb$POS, vcf.comb$REF, vcf.comb$ALT, sep='_')
vcf.comb$CADD_phred <- cadd[vcf.uniq.var, "PHRED"]

# remove rows where all samples are zero
vcf.comb <- vcf.comb[rowSums(gt.agg , na.rm = T) >0, ]

# change col names so cadd is earlier
col.order <- c(1:14, ncol(vcf.comb), 15:(ncol(vcf.comb)-1))
vcf.comb <- vcf.comb[, col.order]

# sort by decreasing CADD_phred
vcf.comb <- vcf.comb[order(vcf.comb$CADD_phred, decreasing=T) ,]

# write out one version with the filter applied, and one without
write.table(vcf.comb, vcf.output.unfilt.f, sep='\t', quote=F, row.names = F, col.names = T)

# keep these GATK filter fileds
keep.gatk.filters <- c('PASS')
write.table(vcf.comb[vcf.comb$FILTER %in% keep.gatk.filters,], vcf.output.filt.f, sep='\t', quote=F, row.names = F, col.names = T)
