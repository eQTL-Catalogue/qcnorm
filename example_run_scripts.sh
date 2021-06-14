
# Runs pop_assign workflow only
nextflow run main.nf -profile tartu_hpc -resume\
 -entry pop_assign_only\
 --study_name CommonMind\
 --vcf_file /gpfs/space/projects/CommonMind/genotypes/imputed/postimpute/crossmap_vcf/CommonMind_genotyped.MAF001.vcf.gz\
 --outdir ./results_test_CommonMind

# Runs QC workflow only
# (pop_assign have been run before and projections.tsv file is provided via pop_assign_projections parameter)
nextflow run main.nf -profile tartu_hpc -resume\
 -entry qc_only\
 --study_name CommonMind\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/CommonMind\
 --pop_assign_projections /gpfs/space/home/kerimov/qcnorm_eqtlcat/results_test_CommonMind/CommonMind/pop_assign/projections_comb.tsv\
 --sample_meta_path /gpfs/space/home/kerimov/SampleArcheology/studies/cleaned/CommonMind.tsv\
 --outdir ./results_test_CommonMind

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