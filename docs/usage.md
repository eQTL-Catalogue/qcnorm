# Normalisation

## Run the normalisation pipeline on a microarray dataset

The `--vcf_file' paramter is only required if you want to use the normalisation output directly in the qtlmap pipeline.

```
nextflow run normalisation.nf -profile tartu_hpc -resume\
 --is_microarray\
 --microarray_exp_matrix_path /gpfs/hpc/projects/eQTLCatalogue/processed/expression_matrices/HumanHT-12_V4/raw/CEDAR.tsv.gz\
 --study_name CEDAR\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/CEDAR.tsv\
 --outdir CEDAR\
 --vcf_file /gpfs/hpc/projects/genomic_references/CEDAR/genotypes/Michigan_GRCh37_1KGPhase3_220918/GRCh38/CEDAR_GRCh38.filtered.renamed.vcf.gz
 ```
