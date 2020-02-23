nextflow run main_proto.nf -profile tartu_hpc\
 --study_name Alasoo_2018\
 --quant_results_path /gpfs/hpc/home/a72094/datasets/processed/Alasoo_2018\
 --sample_meta_path /gpfs/hpc/home/kerimov/SampleArcheology/studies/cleaned/Alasoo_2018_test.tsv\
 --skip_leafcutter_norm\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_txrev_norm

