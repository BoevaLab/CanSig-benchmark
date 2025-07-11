FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC

# Install system dependencies
RUN apt-get update && apt-get install -y \
    # Basic utilities
    wget \
    curl \
    git \
    build-essential \
    software-properties-common \
    apt-transport-https \
    ca-certificates \
    gnupg \
    lsb-release \
    # Development tools
    gcc \
    g++ \
    gfortran \
    make \
    cmake \
    pkg-config \
    # Libraries for R packages
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libssh2-1-dev \
    # Libraries for Python packages
    libffi-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    libncurses5-dev \
    libncursesw5-dev \
    xz-utils \
    tk-dev \
    libgdbm-dev \
    libc6-dev \
    libbz2-dev \
    libexpat1-dev \
    liblzma-dev \
    zlib1g-dev \
    # HDF5 and related
    libhdf5-dev \
    libnetcdf-dev \
    libopenblas-dev \
    liblapack-dev \
    # Geospatial libraries
    libgdal-dev \
    libproj-dev \
    libgeos-dev \
    libudunits2-dev \
    # Graphics libraries
    libcairo2-dev \
    libxt-dev \
    libx11-dev \
    libxext-dev \
    # Additional dependencies
    pandoc \
    pandoc-citeproc \
    && rm -rf /var/lib/apt/lists/*

# Install R 4.3.3 first (before changing Python)
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" && \
    apt-get update && \
    apt-get install -y r-base r-base-dev && \
    rm -rf /var/lib/apt/lists/*

# Install Python 3.11
RUN add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y \
    python3.11 \
    python3.11-dev \
    python3.11-distutils \
    python3.11-venv \
    && rm -rf /var/lib/apt/lists/*

# Install pip for Python 3.11 (don't change system Python defaults)
RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.11

# Create symlinks for convenience but keep system Python intact
RUN ln -sf /usr/bin/python3.11 /usr/local/bin/python && \
    ln -sf /usr/bin/python3.11 /usr/local/bin/python3

# Install Python packages
RUN pip3 install --no-cache-dir \
    numpy==1.24.4 \
    pandas==2.1.0 \
    scipy==1.11.2 \
    scikit-learn==1.3.0 \
    matplotlib==3.7.2 \
    seaborn==0.12.2 \
    anndata==0.9.2 \
    scanpy==1.7.2 \
    umap-learn==0.5.3 \
    numba==0.57.1 \
    h5py==3.11.0 \
    tables==3.9.2 \
    statsmodels==0.14.0 \
    joblib==1.3.2 \
    networkx==3.1 \
    pillow==10.3.0 \
    requests==2.31.0 \
    pyyaml==6.0.1 \
    tqdm==4.66.1 \
    packaging==23.1 \
    patsy==0.5.3 \
    cycler==0.11.0 \
    kiwisolver==1.4.5 \
    pyparsing==3.0.9 \
    python-dateutil==2.8.2 \
    pytz==2023.3.post1 \
    six==1.16.0 \
    threadpoolctl==3.2.0 \
    fonttools==4.42.1 \
    contourpy==1.1.0 \
    munkres==1.0.7 \
    certifi==2024.8.30 \
    charset-normalizer==3.2.0 \
    urllib3==2.0.4 \
    idna==3.4 \
    setuptools==68.1.2 \
    wheel==0.41.2 \
    natsort==8.4.0 \
    pynndescent==0.5.10 \
    llvmlite==0.40.1 \
    numexpr==2.8.4 \
    blosc==1.11.1 \
    lz4 \
    snappy \
    brotli==1.1.0

# Set up R library path
ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN mkdir -p $R_LIBS_USER

# Install BiocManager first
RUN R -e "install.packages('BiocManager', repos='https://cran.rstudio.com/')"

# Install essential R packages first
RUN R -e "install.packages(c('devtools', 'reticulate', 'remotes'), repos='https://cran.rstudio.com/')"

# Install Bioconductor packages
RUN R -e "BiocManager::install(c( \
    'Biobase', \
    'BiocGenerics', \
    'S4Vectors', \
    'IRanges', \
    'GenomeInfoDb', \
    'GenomicRanges', \
    'SummarizedExperiment', \
    'SingleCellExperiment', \
    'DelayedArray', \
    'HDF5Array', \
    'beachmat', \
    'BiocNeighbors', \
    'BiocParallel', \
    'BiocSingular', \
    'scuttle', \
    'batchelor', \
    'rhdf5', \
    'zlibbioc' \
    ), update=FALSE)"

# Install CRAN packages
RUN R -e "install.packages(c( \
    'Matrix', \
    'Rcpp', \
    'RcppArmadillo', \
    'RcppEigen', \
    'ggplot2', \
    'dplyr', \
    'tidyr', \
    'stringr', \
    'magrittr', \
    'tibble', \
    'readr', \
    'purrr', \
    'forcats', \
    'lubridate', \
    'glue', \
    'rlang', \
    'cli', \
    'crayon', \
    'pillar', \
    'vctrs', \
    'lifecycle', \
    'fansi', \
    'utf8', \
    'ellipsis', \
    'digest', \
    'jsonlite', \
    'curl', \
    'mime', \
    'R6', \
    'fastmap', \
    'htmltools', \
    'htmlwidgets', \
    'igraph', \
    'plotly', \
    'RColorBrewer', \
    'viridis', \
    'scales', \
    'gridExtra', \
    'cowplot', \
    'patchwork', \
    'pheatmap', \
    'cluster', \
    'MASS', \
    'survival', \
    'nlme', \
    'mgcv', \
    'lattice', \
    'KernSmooth', \
    'class', \
    'e1071', \
    'randomForest', \
    'rpart', \
    'nnet', \
    'foreign', \
    'boot', \
    'spatial', \
    'snow', \
    'parallel', \
    'logger', \
    'optparse' \
    ), repos='https://cran.rstudio.com/')"

# Install Seurat and dependencies

# Install monocle3 (requires special handling)
RUN R -e "devtools::install_github('cole-trapnell-lab/leidenbase')" && \
    R -e "devtools::install_github('cole-trapnell-lab/monocle3')"

# Configure reticulate to use system Python
RUN R -e "reticulate::use_python('/usr/bin/python3.11', required=TRUE)"

# Set environment variables for reticulate
ENV RETICULATE_PYTHON=/usr/bin/python3.11
ENV PYTHONPATH=/usr/lib/python3.11/site-packages

RUN R -e "remotes::install_version('Seurat', version = '4.3.0', repos = 'https://cran.rstudio.com/')"
RUN R -e "devtools::install_github('cellgeni/sceasy')"
ENV R_MAKEVARS_USER=/tmp/Makevars
RUN echo "CFLAGS += -Wno-error=format-security -Wno-format-extra-args -Wno-format" > /tmp/Makevars && \
    echo "CXXFLAGS += -Wno-error=format-security -Wno-format-extra-args -Wno-format" >> /tmp/Makevars
RUN R -e "devtools::install_github('https://github.com/BoevaLab/scalop.git')"
RUN R -e 'devtools::install_version("SeuratObject", version = "4.1.3", dependencies = TRUE, upgrade = "never")'
RUN R -e 'devtools::install_version("Seurat", version = "4.3.0", dependencies = TRUE, upgrade = "never")'
RUN R -e 'devtools::install_version("GeneNMF", version = "0.4.0", dependencies = TRUE, upgrade = "never")'


# Set working directory
WORKDIR /workspace

# Set default command
CMD ["R"]
