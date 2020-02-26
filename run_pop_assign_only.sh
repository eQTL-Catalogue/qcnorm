nextflow run pop_assign.nf -profile tartu_hpc -resume\
 --study_name Alasoo_2018\
 --vcf_file /gpfs/hpc/home/kerimov/qcnorm/data/Alasoo_2018_filtered_chr21_test_4samples.vcf.gz\
 --ref_genome /gpfs/hpc/home/kerimov/1000G_genome/source_data/GRCh38_renamed_ids_no_multiallelic.vcf.gz\
 --populations_file /gpfs/hpc/home/kerimov/1000G_genome/source_data/1000G_sample_metadata.tsv\
 --exclude_population\
 --ids_to_remove_file /gpfs/hpc/home/kerimov/1000G_genome/source_data/amrs.txt\
 --num_pc 3

# Mandatory parameters: study_name, vcf_file, ref_genome, populations_file
 

