nextflow run main.nf -profile tartu_hpc -resume\
 --study_name Alasoo_2018\
 --vcf_file /gpfs/hpc/home/kerimov/qcnorm/data/Alasoo_2018_filtered_chr21_test_4samples.vcf.gz\
 --exclude_population\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/processed/Alasoo_2018\
 --sample_meta_path /gpfs/hpc/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018_test.tsv\
 --outdir ./results/comp_pipeline_results
