nextflow run quality_control.nf -profile tartu_hpc -resume \
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/processed/Alasoo_2018\
 --study_name Alasoo_2018\
 --sample_meta_path /gpfs/hpc/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018_test.tsv\
 --pop_assign_projections /gpfs/hpc/home/kerimov/qcnorm/results/Alasoo_2018/pop_assign/projections_comb.tsv
#  --ge_exp_matrix_path /gpfs/hpc/projects/eQTLCatalogue/processed/Alasoo_2018/featureCounts/merged_gene_counts.txt\
#  --mbv_files_dir /gpfs/hpc/projects/eQTLCatalogue/processed/Alasoo_2018/MBV\



 