message(" ## Loading libraries: optparse")
suppressPackageStartupMessages(library("optparse"))


#Parse command-line options
option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  make_option(c("-c", "--count_matrix"), type="character", default=NULL,
              help="Counts matrix file path. Tab separated file", metavar = "type"),
  make_option(c("-s", "--sample_meta"), type="character", default=NULL,
              help="Sample metadata file. Tab separated file", metavar = "type"),
  make_option(c("-p", "--phenotype_meta"), type="character", default=NULL,
              help="Phenotype metadata file. Tab separated file", metavar = "type"),
  make_option(c("-q", "--quant_method"), type="character", default="gene_counts",
              help="Quantification method. Possible values: gene_counts, leafcutter, txrevise, transcript_usage, exon_counts and HumanHT-12_V4 [default \"%default\"]", metavar = "type"),
  make_option(c("-o", "--outdir"), type="character", default="./normalised_results/",
              help="Path to the output directory. [default \"%default\"]", metavar = "type"),
  make_option(c("-n", "--name_of_study"), type="character", default=NULL,
              help="Name of the study. Optional", metavar = "type"),
  make_option(c("-t", "--tpm_quantile_file"), type="character", default=NULL,
              help="TPM quantile TSV file with phenotype_id column", metavar = "type"),
  make_option(c("--filter_qc"), type="logical", default=FALSE,
              help="Flag to filter out samples that have failed QC [default \"%default\"]", metavar = "bool"),
  make_option(c("--keep_XY"), type="logical", default=FALSE,
              help="Keep genes on the X and Y chromosomes [default \"%default\"]", metavar = "bool"),
  make_option(c("--eqtlutils"), type="character", default=NULL,
              help="Optional path to the eQTLUtils R package location. If not specified then eQTLUtils is assumed to be installed in the container. [default \"%default\"]", metavar = "type")
)

message(" ## Parsing options")
opt <- optparse::parse_args(OptionParser(option_list=option_list))

message(" ## Loading libraries: devtools, dplyr, SummarizedExperiment, cqn, data.table")
suppressPackageStartupMessages(library("devtools"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("cqn"))
suppressPackageStartupMessages(library("data.table"))

#Debugging
if (FALSE) {
  opt = list()
  opt$n = "GTEx"
  opt$c="TPM/"
  opt$s="GTEx.tsv"
  opt$p="txrevise_Ensembl_96_phenotype_metadata.tsv.gz"
  opt$q="txrevise"
  opt$o="."
  opt$t="GTEx_95quantile_tpm.tsv.gz"
  opt$keep_XY=TRUE
  opt$filter_qc=TRUE
}

count_matrix_path = opt$c
sample_meta_path = opt$s
phenotype_meta_path = opt$p
output_dir = opt$o
quant_method = opt$q
study_name = opt$n
quantile_tpm_path = opt$t
filter_qc = opt$filter_qc
eqtlutils_path = opt$eqtlutils
keep_XY = opt$keep_XY
tpm_threshold = 1

message("######### Options: ######### ")
message("######### Working Directory  : ", getwd())
message("######### quant_method       : ", quant_method)
message("######### count_matrix_path  : ", count_matrix_path)
message("######### sample_meta_path   : ", sample_meta_path)
message("######### phenotype_meta_path: ", phenotype_meta_path)
message("######### output_dir         : ", output_dir)
message("######### opt_study_name     : ", study_name)
message("######### filter_qc          : ", filter_qc)
message("######### eqtlutils_path     : ", eqtlutils_path)
message("######### keep_XY            : ", keep_XY)
message("######### tpm_threshold      : ", tpm_threshold)


#Load eQTLUtils
if (!is.null(eqtlutils_path)){
  devtools::load_all(eqtlutils_path)
}

dummy <- assertthat::assert_that(!is.null(count_matrix_path) && file.exists(count_matrix_path), msg = paste0("count_matrix_path: \"", count_matrix_path, "\" is missing"))
dummy <- assertthat::assert_that(!is.null(sample_meta_path) && file.exists(sample_meta_path), msg = paste0("sample_meta_path: \"", sample_meta_path, "\" is missing"))
dummy <- assertthat::assert_that(!is.null(phenotype_meta_path) && file.exists(phenotype_meta_path), msg =paste0("phenotype_meta_path: \"", phenotype_meta_path, "\" is missing"))

# Read the inputs
message("## Reading sample metadata ##")
sample_metadata <- utils::read.csv(sample_meta_path, sep = '\t', stringsAsFactors = FALSE)

if (is.null(study_name)) { 
  assertthat::has_name(sample_metadata, "study" )
  study_name <- sample_metadata$study[1] 
}

message("## Reading phenotype metadata ##")
phenotype_meta = readr::read_delim(phenotype_meta_path, delim = "\t", col_types = "ccccciiicciidi")
quantile_tpm_df = NULL
if(!is.null(quantile_tpm_path)) {
  quantile_tpm_df = readr::read_delim(quantile_tpm_path, delim = "\t", col_types = "cccd")
}

message("## Reading count matrix ##")
if (quant_method == "txrevise") {
  data_fc <- eQTLUtils::importTxreviseCounts(count_matrix_path)
} else if (quant_method == "leafcutter") {
  #Use any number of white spaces for leafcutter data (mixed spaces and tabs in the header, needs to be fixed)
  data_fc <- utils::read.csv(count_matrix_path, sep = "", stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data_fc <- utils::read.csv(count_matrix_path, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
}

message("## Make Summarized Experiment ##")
se <- eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = data_fc, row_data = phenotype_meta, col_data = sample_metadata, quant_method = quant_method)
dim(se)

if (filter_qc){
  message("## Filter SummarizedExperiment by removing samples that fail QC ##")
  se <- eQTLUtils::filterSummarizedExperiment(se, filter_rna_qc = TRUE, filter_genotype_qc = TRUE)
  dim(se)
}

if (!dir.exists(paste0(output_dir, "/qtl_group_split_norm"))){
  dir.create(paste0(output_dir, "/qtl_group_split_norm"), recursive = TRUE)
}
if (!dir.exists(paste0(output_dir, "/norm_not_filtered"))){
  dir.create(paste0(output_dir, "/norm_not_filtered"), recursive = TRUE)
}

split_and_filter_by_qtlgroup <- function(norm_count_df, 
                                         sample_metadata_df, 
                                         phenotype_metadata,
                                         qtl_group_selected,
                                         study_name,
                                         out_dir, 
                                         quantile_tpms, 
                                         tpm_thres) {
  message("Processing qtl_group: ", qtl_group_selected)
  
  # extract genes which are expressed (TPM > tpm threshold)
  expressed_genes <- quantile_tpms %>% 
    dplyr::filter(qtl_group == qtl_group_selected, median_tpm > tpm_thres) %>% 
    dplyr::pull(phenotype_id)
  message("Expressed gene count by quantile TPM: ", length(expressed_genes))
  
  # extract the phenotype metadata of expressed genes
  expressed_phenotype_metadata <- phenotype_metadata %>% 
    dplyr::filter(gene_id %in% expressed_genes)
  
  # extract sample_id's belonging to selected qtl_group_selected
  qtl_group_samples = sample_metadata_df %>% 
    dplyr::as_tibble() %>% 
    dplyr::filter(qtl_group == qtl_group_selected) %>% 
    dplyr::pull(sample_id)
  
  message("Samples in qtl_group: ", length(qtl_group_samples))
  
  # extract qtl_group samples and expressed phenotypes from normalised count matrix
  counts_for_selected_qtl_group <- norm_count_df %>% 
    dplyr::as_tibble() %>% 
    dplyr::select(c("phenotype_id", qtl_group_samples)) %>% 
    dplyr::filter(phenotype_id %in% expressed_phenotype_metadata$phenotype_id)
  
  message("Count of phenotypes for ", qtl_group_selected, " qtl_group is: ", nrow(counts_for_selected_qtl_group))
  message("Count of samples for ", qtl_group_selected, " qtl_group is: ", ncol(counts_for_selected_qtl_group) - 1)
  utils::write.table(x = counts_for_selected_qtl_group, 
                     file = file.path(out_dir, "qtl_group_split_norm", paste0(study_name, ".", qtl_group_selected ,".tsv")), 
                     sep = "\t", quote = FALSE, 
                     row.names = FALSE, col.names = TRUE)
  message("TPM filtered and splitted by qtl_group TSVs exported into: ", file.path(out_dir, "qtl_group_split_norm", paste0(study_name, ".", qtl_group_selected ,".tsv")))
  message("______________________________________________________")
}

# Get unique qtl_group values
qtl_groups_in_se <- se$qtl_group %>% base::unique()

# Show sample sizes of each qtl_group
table(se$qtl_group)

# Create empty dataframes to collect quantile and median TPM values from 
# different qtl_groups within the study
merged_median_tpm_df <- base::data.frame()
merged_95quantile_tpm_df <- base::data.frame()

for (qtl_group_in_se in qtl_groups_in_se) {
  message("Normalising qtl_group: ", qtl_group_in_se)
  
  # Extract the se of selected qtl_group (qtl_group_in_se)
  se_qtl_group <- se[,se$qtl_group == qtl_group_in_se]
  dim(se_qtl_group)
  
  message("## Starting normalisation process... ##")
  if (quant_method=="gene_counts") {
    cqn_norm <- eQTLUtils::qtltoolsPrepareSE(se_qtl_group, "gene_counts", filter_genotype_qc = FALSE, filter_rna_qc = FALSE, keep_XY)
    cqn_int_norm <- eQTLUtils::normaliseSE_quantile(cqn_norm, assay_name = "cqn")
    
    cqn_int_assay_fc_formatted <- SummarizedExperiment::cbind(phenotype_id = rownames(assays(cqn_int_norm)[["qnorm"]]), assays(cqn_int_norm)[["qnorm"]])
    utils::write.table(cqn_int_assay_fc_formatted, file.path(output_dir, "norm_not_filtered", paste0(study_name, ".", qtl_group_in_se, ".gene_counts_cqn_int_norm.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    message("## Normalised gene count matrix exported into: ", file.path(output_dir, "norm_not_filtered", paste0(study_name, ".", qtl_group_in_se, ".gene_counts_cqn_int_norm.tsv")))
    
    message("## Caclulate median TPM in each biological context ##")
    median_tpm_df = eQTLUtils::estimateMedianTPM(cqn_norm, subset_by = "qtl_group", assay_name = "cqn", prob = 0.5)
    # gzfile = gzfile(file.path(output_dir, paste0(study_name ,"_median_tpm.tsv.gz")), "w")
    # write.table(median_tpm_df, gzfile, sep = "\t", row.names = F, quote = F)
    # close(gzfile)
    merged_median_tpm_df <- merged_median_tpm_df %>% base::rbind(median_tpm_df)
    
    message("## Caclulate 95% quantile TPM in each biological context ##")
    quantile_tpm_df = eQTLUtils::estimateMedianTPM(cqn_norm, subset_by = "qtl_group", assay_name = "cqn", prob = 0.95)
    # gzfile = gzfile(file.path(output_dir, paste0(study_name ,"_95quantile_tpm.tsv.gz")), "w")
    # write.table(quantile_tpm_df, gzfile, sep = "\t", row.names = F, quote = F)
    # close(gzfile)
    # message("## Median tpm values matrix exported into: ", file.path(output_dir, paste0(study_name ,"_median_tpm.tsv.gz")))
    merged_95quantile_tpm_df <- merged_95quantile_tpm_df %>% base::rbind(quantile_tpm_df)
    
    split_and_filter_by_qtlgroup(norm_count_df = cqn_int_assay_fc_formatted, 
                                 sample_metadata = SummarizedExperiment::colData(cqn_int_norm), 
                                 phenotype_metadata = phenotype_meta, 
                                 qtl_group_selected = qtl_group_in_se,
                                 study_name = study_name, 
                                 out_dir = output_dir, 
                                 quantile_tpms = quantile_tpm_df, 
                                 tpm_thres = tpm_threshold)
  } else if (quant_method=="exon_counts") {
    cqn_norm <- eQTLUtils::qtltoolsPrepareSE(se_qtl_group, "exon_counts", filter_genotype_qc = FALSE, filter_rna_qc = FALSE, keep_XY)
    cqn_int_norm <- eQTLUtils::normaliseSE_quantile(cqn_norm, assay_name = "cqn")
    cqn_int_assay_fc_formatted <- SummarizedExperiment::cbind(phenotype_id = rownames(assays(cqn_int_norm)[["qnorm"]]), assays(cqn_int_norm)[["qnorm"]])
    utils::write.table(cqn_int_assay_fc_formatted, file.path(output_dir, "norm_not_filtered", paste0(study_name, ".", qtl_group_in_se ,".exon_counts_cqn_int_norm.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    message("## Normalised exon count matrix exported into: ", output_dir, "/norm_not_filtered/", study_name, ".", qtl_group_in_se ,".exon_counts_cqn_int_norm.tsv")
    
    split_and_filter_by_qtlgroup(norm_count_df = cqn_int_assay_fc_formatted, 
                                 sample_metadata = SummarizedExperiment::colData(cqn_int_norm), 
                                 phenotype_metadata = phenotype_meta,
                                 qtl_group_selected = qtl_group_in_se,
                                 study_name = study_name, 
                                 out_dir = output_dir, 
                                 quantile_tpms = quantile_tpm_df, 
                                 tpm_thres = tpm_threshold)
  } else if (quant_method %in% c("transcript_usage", "txrevise")) {
    q_norm <- eQTLUtils::qtltoolsPrepareSE(se_qtl_group, "txrevise", filter_genotype_qc = FALSE, filter_rna_qc = FALSE, keep_XY)
    qnorm_assay_fc_formatted <- SummarizedExperiment::cbind(phenotype_id = rownames(assays(q_norm)[["qnorm"]]), assays(q_norm)[["qnorm"]])
    utils::write.table(qnorm_assay_fc_formatted, file.path(output_dir, "norm_not_filtered", paste0(study_name, ".", qtl_group_in_se , "." , quant_method, "_qnorm.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)    
    message("## Normalised transcript usage matrix exported into: ", file.path(output_dir, "norm_not_filtered", paste0(study_name, ".", qtl_group_in_se , "." , quant_method, "_qnorm.tsv")))
    
    split_and_filter_by_qtlgroup(norm_count_df = qnorm_assay_fc_formatted, 
                                 sample_metadata = SummarizedExperiment::colData(q_norm), 
                                 phenotype_metadata = phenotype_meta,
                                 qtl_group_selected = qtl_group_in_se,
                                 study_name = study_name, 
                                 out_dir = output_dir, 
                                 quantile_tpms = quantile_tpm_df, 
                                 tpm_thres = tpm_threshold)
  } else if (quant_method == "leafcutter") {
    q_norm <- eQTLUtils::qtltoolsPrepareSE(se_qtl_group, "leafcutter", filter_genotype_qc = FALSE, filter_rna_qc = FALSE, keep_XY)
    qnorm_assay_fc_formatted <- SummarizedExperiment::cbind(phenotype_id = rownames(assays(q_norm)[["qnorm"]]), assays(q_norm)[["qnorm"]])
    utils::write.table(qnorm_assay_fc_formatted, file.path(output_dir, "norm_not_filtered", paste0(study_name, ".", qtl_group_in_se , "." , quant_method, "_qnorm.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)    
    message("## Normalised transcript usage matrix exported into: ", file.path(output_dir, "norm_not_filtered", paste0(study_name, ".", qtl_group_in_se , "." , quant_method, "_qnorm.tsv")))
    
    split_and_filter_by_qtlgroup(norm_count_df = qnorm_assay_fc_formatted, 
                                 sample_metadata = SummarizedExperiment::colData(q_norm), 
                                 phenotype_metadata = phenotype_meta,
                                 qtl_group_selected = qtl_group_in_se,
                                 study_name = study_name, 
                                 out_dir = output_dir, 
                                 quantile_tpms = quantile_tpm_df, 
                                 tpm_thres = tpm_threshold)
  } else if (quant_method == "HumanHT-12_V4") {
    q_norm <- eQTLUtils::qtltoolsPrepareSE(se_qtl_group, "HumanHT-12_V4", filter_genotype_qc = FALSE, filter_rna_qc = FALSE, keep_XY)
    q_int_norm <- eQTLUtils::normaliseSE_quantile(q_norm, assay_name = "norm_exprs")
    qnorm_assay_fc_formatted <- SummarizedExperiment::cbind(phenotype_id = rownames(assays(q_int_norm)[["qnorm"]]), assays(q_int_norm)[["qnorm"]])
    utils::write.table(qnorm_assay_fc_formatted, file.path(output_dir, "qtl_group_split_norm", paste0(study_name, ".", qtl_group_in_se , ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    message("## Normalised HumanHT-12_V4 matrix exported to: ", file.path(output_dir, "qtl_group_split_norm", paste0(study_name, ".", qtl_group_in_se , ".tsv")))
  }
}

if (quant_method=="gene_counts") {
  gzfile = gzfile(file.path(output_dir, paste0(study_name ,"_median_tpm.tsv.gz")), "w")
  write.table(merged_median_tpm_df, gzfile, sep = "\t", row.names = F, quote = F)
  close(gzfile)
  message("## Median TPM values matrix exported into: ", file.path(output_dir, paste0(study_name ,"_median_tpm.tsv.gz")))
  
  gzfile = gzfile(file.path(output_dir, paste0(study_name ,"_95quantile_tpm.tsv.gz")), "w")
  write.table(merged_95quantile_tpm_df, gzfile, sep = "\t", row.names = F, quote = F)
  close(gzfile)
  message("## Quantile TPM values matrix exported into: ", file.path(output_dir, paste0(study_name ,"_95quantile_tpm.tsv.gz")))
}
