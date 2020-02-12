nextflow run main_proto.nf -profile tartu_hpc,pop_assign\
 --is_microarray\
 --exp_matrix_path /gpfs/hpc/home/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/raw/Naranbhai_2015.tsv.gz\
 --study_name Naranbhai_2015\
 --sample_meta_path /gpfs/hpc/home/kerimov/SampleArcheology/studies/cleaned/Naranbhai_2015.tsv

 