# Package configuration for May 2022 to ensure stack compatibility (Ubuntu 20.04, Python 3.8 and R 4.2.0) 
FROM ubuntu:20.04

# Noninteractive since we're running in Docker
ENV DEBIAN_FRONTEND=noninteractive

# Install general dependencies
RUN apt-get update -y && apt-get install -y \
    build-essential gfortran make cmake \
    libcurl4-openssl-dev libxml2-dev libssl-dev \
    curl \
    libopenblas-dev liblapack-dev \
    dirmngr gnupg \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# OpenBLAS makes some of the python and R computation more efficient (particularly with ML)
ENV OPENBLAS_NUM_THREADS=4 

# Install python 3.8 specifically
RUN apt-get update -y && apt-get install -y \
    python3.8 python3.8-venv \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Set up python and any required packages
RUN python3.8 -m venv miRsight-venv
ENV PATH /miRsight-venv/bin:$PATH
RUN . activate
RUN pip3.8 install scikit-learn==1.1.1 pandas==1.4.2

# Install R 4.2.0 specifically
RUN curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/cran.gpg] https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/" | tee /etc/apt/sources.list.d/cran.list \
    && apt-get update 
RUN apt-get install -y r-base=4.2.0-* r-base-dev=4.2.0-* r-base-html=4.2.0-* r-recommended=4.2.0-* \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up R and any required packages
RUN mkdir -p /rlib && chown -R root:root /rlib
ENV R_LIBS_USER=/rlib
RUN R -e "install.packages(c('jsonlite', 'stringi', 'stringr', 'dplyr', 'mice'), lib='/rlib')"
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', lib='/rlib')" && \
    R -e "BiocManager::install(c('ensembldb', 'AnnotationHub', 'BSgenome.Hsapiens.NCBI.GRCh38', 'GenomicScores'), lib='/rlib')"

# Install the ViennaRNA Suite 2.5.0 specifically
RUN curl -L -O https://github.com/ViennaRNA/ViennaRNA/releases/download/v2.5.0/ViennaRNA-2.5.0.tar.gz
RUN tar -zxvf ViennaRNA-2.5.0.tar.gz
RUN rm -f ViennaRNA-2.5.0.tar.gz
WORKDIR /ViennaRNA-2.5.0
RUN ./configure
RUN make
RUN make install

# Copy over app code
WORKDIR /app
COPY . .

# Run
CMD ["python3.8", "-u", "main.py"]