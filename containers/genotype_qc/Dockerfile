FROM nfcore/base
LABEL authors="kaur.alasoo@ut.ee" \
      description="Docker image containing all requirements for the kauralasoo/genotype_qc pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-rnaseq-1.3/bin:$PATH
COPY bin/ldak /usr/bin/
COPY bin/GenotypeHarmonizer-1.4.20/ /usr/bin/
