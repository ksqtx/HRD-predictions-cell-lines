FROM rocker/tidyverse:4.1.3

RUN apt-get update && apt-get install -y \
  ## For Miniconda
  bzip2 \
  wget \
  ## For Rhtslib package
  libbz2-dev \
  liblzma-dev

# Install Miniconda
ENV CONDA_DIR=/opt/conda

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
  && /bin/bash Miniconda3-latest-Linux-x86_64.sh -bfp ${CONDA_DIR} \
  && rm -rf Miniconda3-latest-Linux-x86_64.sh

## add conda to PATH
ENV PATH=${CONDA_DIR}/bin:$PATH

## add channels
RUN conda config --add channels bioconda \
  && conda config --add channels conda-forge

# Install bcftools
RUN conda update -y conda \
  && conda install -y \
    -c bioconda \
    bcftools=1.15.1

# Install R packages
RUN R -e "install.packages(c('devtools', 'randomForest', 'scales', 'ranger', 'Boruta', 'doMC', 'BiocManager'))"
RUN R -e "BiocManager::install(pkgs = c('BSgenome', 'BSgenome.Hsapiens.UCSC.hg38', 'StructuralVariantAnnotation', 'VariantAnnotation'), ask = FALSE)"
RUN R -e "devtools::install_github('linxihui/NNLM@4574bca', dependencies = TRUE)"
RUN R -e "devtools::install_github('UMCUGenetics/mutSigExtractor@f218803', dependencies = TRUE)"
RUN R -e "devtools::install_github('UMCUGenetics/CHORD@594916d', dependencies = FALSE)"
