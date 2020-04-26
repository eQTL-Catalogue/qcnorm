# Normalisation

## Run the normalisation pipeline on a microarray dataset

```
nextflow run normalisation.nf -profile tartu_hpc -resume\
 --is_microarray\
 --microarray_exp_matrix_path /gpfs/hpc/projects/eQTLCatalogue/processed/expression_matrices/HumanHT-12_V4/raw/CEDAR.tsv.gz\
 --study_name CEDAR\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/CEDAR.tsv\
 --outdir test_results
 ```
