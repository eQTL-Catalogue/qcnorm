FROM bioconductor/bioconductor_docker:RELEASE_3_15
LABEL authors="Kaur Alasoo" \
      description="Docker image containing all requirements for making sample relatedness heatmap"

RUN R -e "BiocManager::install(c('dplyr', 'optparse', 'readr', 'tidyr','pheatmap'))"
