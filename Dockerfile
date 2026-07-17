FROM rocker/tidyverse:4.4

RUN install2.r \
rmarkdown \
DT \
ggpubr \
pheatmap \
pvclust \
circlize \
littler \
gplots \
here \
ggpmisc \
gt \
Hmisc \
BiocManager \
UpSetR \
dendextend \
hrbrthemes \
optparse \
readr \
ggthemes \
Seurat \
SeuratObject \
Matrix \
harmony \
tibble \
ggtext \
svglite \
systemfonts \
ggforce

RUN sudo apt-get update -y && sudo apt-get install -y \
libglpk-dev libbz2-dev libproj-dev libgdal-dev

RUN R -e "install.packages(c('proj4', 'ggalt'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

# install pak
RUN R -e "install.packages('pak')"

############################################################
# MINICONDA
############################################################
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
chmod +x miniconda.sh && \
./miniconda.sh -b -p /home/rstudio/.local/share/r-miniconda && \
rm miniconda.sh

ENV CONDA_DIR="/home/rstudio/.local/share/r-miniconda"
ENV PATH="${CONDA_DIR}/bin:${PATH}"
ENV CONDA_PLUGINS_AUTO_ACCEPT_TOS=true

# avoid mamba solver completely
RUN conda config --set solver classic

############################################################
# pre-create omnideconv env (no runtime creation)
############################################################
RUN conda create -y -n r-omnideconv python=3.8 pip

# install python dependency that was failing
RUN conda run -n r-omnideconv pip install --no-cache-dir anndata

############################################################
# install omnideconv AFTER env exists
############################################################
RUN R -e "pak::pkg_install('omnideconv/omnideconv')"

RUN Rscript -e "BiocManager::install(c('SimBu', 'DESeq2', 'biomaRt', 'ComplexHeatmap'))"

############################################################
# prevent reticulate from trying to manage conda at runtime
############################################################
ENV RETICULATE_MINICONDA_ENABLED=FALSE
ENV RETICULATE_PYTHON=${CONDA_DIR}/envs/r-omnideconv/bin/python