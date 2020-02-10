#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
options(bitmapType='cairo')

library(ggplot2)
library(dplyr)
library(DMwR)

pca_table = if (is.na(args[1])) 'main_overlapped_pca.vect' else args[1]
projections = if (is.na(args[2])) 'new_dataset_scores.profile.adj' else args[2]
source_populations_file = if (is.na(args[3])) "samples_data.tsv" else args[3]
data_name = if (is.na(args[4])) 'New_dataset' else args[4] 
output_dir = if (is.na(args[5])) 'plots' else args[5] 

save_ggplots <- function(plot, path = "plots", filename = "unnamed_plot", height = 15, width = 25){
  if (!dir.exists(path)){
    dir.create(path, recursive = TRUE)
  }
  ggsave(plot = plot,
         filename = paste0(filename, ".png"), 
         path = path,
         device = "png", 
         height = height, 
         width = width,
         units = "cm",
         dpi = 300)
  
  ggsave(plot = plot,
         filename = paste0(filename, ".pdf"), 
         path = path,
         device = "pdf", 
         height = height, 
         width = width,
         units = "cm",
         dpi = 300)
}

pca <- read.table(pca_table, header = FALSE, )
pca <- pca %>% select(-c(2)) %>% rename(genotype_id=V1, PC1=V3, PC2=V4, PC3=V5)

source_populations = read.table(source_populations_file, header = TRUE, sep='\t')
main_pca <- merge(pca, select(source_populations, genotype_id, superpopulation_code), by  = "genotype_id") 

ref_pca_plot <- ggplot(main_pca, aes(x=PC1, y=PC2, color=as.factor(superpopulation_code))) +
  geom_point() + guides(color = guide_legend(title='Superpopulation')) + 
  labs(x = "PC1", y = "PC2") + ggtitle(data_name)+ coord_fixed()

save_ggplots(plot = ref_pca_plot, filename = "ref_pca")

projections_pcs <- read.table(projections, header = TRUE)
projections_pcs <- projections_pcs %>% select(-c(2,3,4)) %>% 
  rename(genotype_id=ID1, PC1=Adjusted1, PC2=Adjusted2, PC3=Adjusted3) %>% 
  mutate(superpopulation_code = data_name)

projections_only_plot <- ggplot(projections_pcs, aes(x=PC1, y=PC2)) + geom_point()+ coord_fixed()
save_ggplots(plot = projections_only_plot, filename = "projections_only")

comb_pcs = rbind(main_pca, projections_pcs)

projected_plot <- ggplot(comb_pcs, aes(x=PC1, y=PC2, color=superpopulation_code)) + geom_point() + 
  guides(color = guide_legend(title='Superpopulation')) + ggtitle(data_name) +
  ggplot2::theme_bw() 
save_ggplots(plot = projected_plot, filename = "projections_on_ref")
write.table(comb_pcs, file='projections_comb.tsv', quote = FALSE, row.names = FALSE, sep = '\t')
################### Populations assignment ################### 
means_main_pca <- main_pca %>% group_by(superpopulation_code) %>% summarize(mean_P1 = mean(PC1), mean_P2 = mean(PC2))

distance <- function(p1, p2, mean1, mean2){
  return(sqrt((mean1-p1)**2 + (mean2-p2)**2))
}

populations = unique(main_pca$superpopulation_code)
samples_with_dists = projections_pcs

for (population_code in populations) {
  col_name = sprintf('dist_%s', population_code)
  samples_with_dists = samples_with_dists %>% 
    mutate( !!col_name := distance(PC1, PC2, means_main_pca[means_main_pca$superpopulation_code == population_code, ]$mean_P1,
                                   means_main_pca[means_main_pca$superpopulation_code == population_code, ]$mean_P2))
}

# samples_with_dists %>% mutate(sum = select(., cols) %>% apply(1, sum))
samples_with_dists = samples_with_dists %>% mutate(sum_dist = select(., starts_with('dist_')) %>% apply(1, sum))

for (population_code in populations) {
  col_name = sprintf('dist_%s', population_code)
  col_name_norm = sprintf('dist_%s_norm', population_code)
  samples_with_dists = samples_with_dists %>%
    mutate( !!col_name_norm := !!as.name(col_name) / sum_dist)
}
############### k-nearest neighbours #############

train <- main_pca %>% select(PC1, PC2, superpopulation_code)
test <- samples_with_dists %>% select(PC1, PC2)

knn_pops <- kNN(superpopulation_code ~ PC1+PC2,train,test,norm=FALSE,k=5)
new_populations <- samples_with_dists %>% mutate(knn_pop = knn_pops)

main_pca_to_map <- main_pca %>% mutate(knn_pop = superpopulation_code)

knn_plot <- ggplot(new_populations, aes(x=PC1, y=PC2, color=as.factor(knn_pop))) + 
  geom_point() + 
  geom_point(data=main_pca_to_map, alpha = 0.3, shape=4) + 
  guides(color = guide_legend(title='KNN population')) + ggtitle(data_name) + coord_fixed() 

save_ggplots(plot = knn_plot, filename = "knn")
threshold = 0.1

new_populations_threshold = new_populations %>% 
  mutate(min_val = apply(select(new_populations, ends_with('norm')), 1, min)) %>% 
  mutate(knn_pop = if_else(min_val>threshold, 'АdMixed', as.character(knn_pop)))  %>% 
  select(-min_val)

knn_threshold <- ggplot(new_populations_threshold, aes(x=PC1, y=PC2, color=as.factor(knn_pop))) + 
  geom_point() + 
  geom_point(data=main_pca_to_map, alpha = 0.3, shape=4) + 
  guides(color = guide_legend(title='KNN population')) + ggtitle(data_name) + coord_fixed() 

save_ggplots(plot = knn_plot, filename = "knn_threshold")

###################################################
finalData <- new_populations %>% mutate(knn_pop_threshold = new_populations_threshold$knn_pop)
write.table(finalData, file='populations.tsv', quote = FALSE, row.names = FALSE, sep = '\t')
