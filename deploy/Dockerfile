FROM bioconductor/bioconductor_docker:RELEASE_3_19

WORKDIR /code

RUN apt-get update && apt-get install -y --no-install-recommends \
    libsqlite3-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    python3-pip \
    python3-venv \
    && apt-get clean

RUN R -e "install.packages(c('shiny','dplyr','ggplot2','tidyr','DT','plotly','RColorBrewer','data.table','DBI','RSQLite','uwot','bslib','circlize'), repos='https://cloud.r-project.org')"

RUN R -e "BiocManager::install(c('clusterProfiler','org.Hs.eg.db'), ask=FALSE, update=FALSE)"

RUN pip install huggingface_hub

COPY app.R .
COPY start.sh .

RUN chmod +x start.sh

CMD ["bash", "start.sh"]