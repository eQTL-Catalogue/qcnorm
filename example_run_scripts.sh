
# Runs pop_assign workflow only
nextflow run main.nf -profile tartu_hpc -resume\
 -entry pop_assign_only\
 --study_name CommonMind\
 --vcf_file /gpfs/space/projects/CommonMind/genotypes/imputed/postimpute/crossmap_vcf/CommonMind_genotyped.MAF001.vcf.gz\
 --outdir ./results_test_CommonMind

# Runs QC workflow only
# (pop_assign have been run before and projections.tsv file is provided via pop_assign_projections parameter)
# (MBV folder under "quant_results_path" is expected to be there with MBV files inside)
nextflow run main.nf -profile tartu_hpc -resume\
 -entry qc_only\
 --study_name CommonMind\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/CommonMind\
 --pop_assign_projections /gpfs/space/home/kerimov/qcnorm_eqtlcat/results_test_CommonMind/CommonMind/pop_assign/projections_comb.tsv\
 --sample_meta_path /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/CommonMind.tsv\
 --outdir ./results_test_CommonMind

# or with a specific path to gene_expression feature counts file
# you can skip the MBV files directory when you don't use "quant_results_path" and use "ge_exp_matrix_path"
nextflow run main.nf -profile tartu_hpc -resume\
 -entry qc_only\
 --study_name LCL_multi\
 --ge_exp_matrix_path /gpfs/space/home/kerimov/multistudy_nf/results/metadata_and_fc_merged/joined_feature_counts.tsv\
 --pop_assign_projections /gpfs/space/home/kerimov/qcnorm/results_multi_LCL_noAMR/LCL_multi/pop_assign/projections_comb.tsv\
 --sample_meta_path /gpfs/space/home/kerimov/multistudy_nf/results/metadata_and_fc_merged/merged_metadata.tsv\
 --outdir ./results_multi_LCL_qc

# Runs normalisation workflow only
nextflow run main.nf -profile tartu_hpc -resume\
 -entry norm_only\
 --study_name CommonMind\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/CommonMind\
 --sample_meta_path /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/CommonMind.tsv\
 --skip_leafcutter_norm\
 --outdir ./results_test_CommonMind

# Runs pop_assign and QC (without normalisation)
nextflow run main.nf -profile tartu_hpc -resume\
 -entry pop_assign_and_qc\
 --study_name CommonMind\
 --vcf_file /gpfs/space/projects/CommonMind/genotypes/imputed/postimpute/crossmap_vcf/CommonMind_genotyped.MAF001.vcf.gz\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/CommonMind\
 --sample_meta_path /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/CommonMind.tsv\
 --outdir ./results_test_CommonMind

# Runs QC and normalisation 
# (pop_assign have been run before and projections.tsv file is provided via pop_assign_projections parameter)
nextflow run main.nf -profile tartu_hpc -resume\
 -entry qc_and_norm\
 --vcf_file /gpfs/space/projects/CommonMind/genotypes/imputed/postimpute/crossmap_vcf/CommonMind_genotyped.MAF001.vcf.gz\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/CommonMind\
 --pop_assign_projections /gpfs/space/home/kerimov/qcnorm_eqtlcat/results_test_CommonMind/CommonMind/pop_assign/projections_comb.tsv\
 --study_name CommonMind\
 --sample_meta_path /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/CommonMind.tsv\
 --skip_leafcutter_norm\
 --outdir ./results_test_CommonMind

# Runs all three workflows: pop_assign, QC and normalisation
nextflow run main.nf -profile tartu_hpc -resume\
 --study_name CommonMind\
 --vcf_file /gpfs/space/projects/CommonMind/genotypes/imputed/postimpute/crossmap_vcf/CommonMind_genotyped.MAF001.vcf.gz\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/CommonMind\
 --sample_meta_path /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/CommonMind.tsv\
 --skip_leafcutter_norm\
 --outdir ./results_test_CommonMind

# Runs normalisation with input.tsv rnaseq
nextflow run main.nf -profile tartu_hpc -resume\
 -entry norm_only_with_tsv\
 --input_tsv /gpfs/space/home/kerimov/qcnorm_fast/input_tsv.tsv\
 --exclude_population\
 --skip_leafcutter_norm\
 --norm_keep_XY\
 --outdir /gpfs/space/home/kerimov/qcnorm_fast/results_GTEx_all
 
 # Runs normalisation for microarray data using tsv input
 nextflow run main.nf -profile tartu_hpc -resume\
   -entry norm_only_with_tsv\
   --input_tsv /gpfs/space/home/kerimov/qcnorm_fast/input_tsv_ma.tsv\
   --is_microarray\
   --norm_keep_XY\
   --outdir eQTL_Catalogue_r5_ma

# Runs normalisation of only RNAseq with a specific count matrix, instead of quant_results_path
# vcf_file is optional (can be skipped, but then qtlmap_inputs table will put "false" for vcf)
nextflow run main.nf -profile tartu_hpc -resume \
  -entry norm_only \
  --study_name GEUVADIS_GBR20 \
  --vcf_file /gpfs/space/projects/eQTLCatalogue/test_data/GEUVADIS_GBR_cohort/vcf/GEUVADIS_GBR20.vcf.gz \
  --ge_exp_matrix_path /gpfs/space/home/kerimov/GitHub/rnaseq/results/test_GEUVADIS_20samples_changed/featureCounts/merged_gene_counts.tsv.gz \
  --sample_meta_path /gpfs/space/projects/eQTLCatalogue/SampleArcheology_V6.0/studies/cleaned/GEUVADIS_test_20samples.tsv \
  --skip_exon_norm \
  --skip_tx_norm \
  --skip_txrev_norm \
  --skip_leafcutter_norm \
  --outdir results/GEUVADIS_GBR20_exp_matrix_norm_only