suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("--kinship"), type="character", default=NULL,
                        help="Path to kinship matrix calculated by gemma.", metavar = "type"),
  optparse::make_option(c("--fam"), type="character", default=NULL,
                        help="Path to PLINK fam file.", metavar = "type"),
  optparse::make_option(c("--out"), type="character", default=NULL,
                        help="Path to the output file.", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Deubgging
if(FALSE){
  opt = list(
    kinship = "testout/GEUVADIS_test_ge.relatedness.sXX.txt",
    fam = "testout/GEUVADIS_test_ge_pruned.fam",
    out = "testout/GEUVADIS_test_ge.kinship.tsv")
}

#Add row and column names to the kinship matrix
kinship_matrix = read.table(opt$kinship)
fam_file = read.table(opt$fam)
colnames(kinship_matrix) = fam_file$V1
kinship_df = cbind(iid = fam_file$V1, kinship_matrix)
write.table(kinship_df, opt$out, sep = "\t", row.names = F, quote = F)
