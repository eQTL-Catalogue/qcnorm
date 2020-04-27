# Popassign

```
nextflow run pop_assign.nf -profile tartu_hpc -resume\
 --study_name GEUVADIS_EUR\
 --vcf_file /gpfs/hpc/projects/genomic_references/GEUVADIS/genotypes/GEUVADIS_GRCh38_filtered.vcf.gz\
 --ref_genome /gpfs/hpc/projects/genomic_references/1000G/pop_assign/GRCh38_renamed_ids_no_multiallelic.vcf.gz\
 --populations_file /gpfs/hpc/projects/genomic_references/1000G/pop_assign/1000G_sample_metadata.tsv\
 --exclude_population\
 --ids_to_remove_file /gpfs/hpc/projects/genomic_references/1000G/pop_assign/amrs.txt\
 --num_pc 3\
 --outdir GEUVADIS_EUR
 ```

# Quality Control

## RNA-seq dataset

```
nextflow run quality_control.nf -profile tartu_hpc -resume \
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/processed/Alasoo_2018\
 --study_name Alasoo_2018\
 --sample_meta_path /gpfs/hpc/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018_test.tsv\
 --pop_assign_projections /gpfs/hpc/home/kerimov/qcnorm/results/Alasoo_2018/pop_assign/projections_comb.tsv
```

# Normalisation

## Mmicroarray dataset (HumanHT-12_V4)

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

## RNA-seq dataset

In this example, we only run gene count and txrevise normalisation, skipping the other three quantification methods. 

```
nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name GEUVADIS_EUR\
 --quant_results_path /gpfs/hpc/home/a72094/projects/rnaseq/results/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/GEUVADIS_EUR.tsv\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --outdir GEUVADIS_EUR
 ```