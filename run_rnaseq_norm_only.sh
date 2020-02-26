nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name Alasoo_2018\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/processed/Alasoo_2018\
 --sample_meta_path /gpfs/hpc/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018_test.tsv\
 --external_quantile_tpm /gpfs/hpc/home/kerimov/qcnorm/results/Alasoo_2018/normalised/Alasoo_2018_95quantile_tpm.tsv.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_txrev_norm\
 --skip_leafcutter_norm
 

