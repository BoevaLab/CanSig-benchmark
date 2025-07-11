# Use Ubuntu 22.04 as base for newer GLIBCXX support
FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies and update libstdc++
RUN apt-get update && apt-get install -y \
    git \
    wget \
    curl \
    build-essential \
    software-properties-common \
    && add-apt-repository ppa:ubuntu-toolchain-r/test \
    && apt-get update \
    && apt-get install -y \
    gcc-11 \
    g++-11 \
    libstdc++6 \
    bzip2 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install Mambaforge
RUN wget -O /tmp/mambaforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh" && \
    bash /tmp/mambaforge.sh -b -p /opt/conda && \
    rm /tmp/mambaforge.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Set working directory
WORKDIR /app

# Copy the environment file
COPY envs/with_build/CanSig-python.yml /tmp/environment.yml

# Create the conda environment
RUN CONDA_OVERRIDE_CUDA="11.8" mamba env create -f /tmp/environment.yml && \
    mamba clean -a --yes

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "CanSig-python", "/bin/bash", "-c"]

# Ensure the environment is activated for subsequent commands
ENV CONDA_DEFAULT_ENV=CanSig-python
ENV PATH=/opt/conda/envs/CanSig-python/bin:$PATH

# Set the default command to python, allowing command line arguments
ENTRYPOINT ["python"]

# Default command (can be overridden)
CMD ["--help"]