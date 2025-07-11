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
    python3.11 \
    python3.11-dev \
    python3.11-distutils \
    && rm -rf /var/lib/apt/lists/*

# Install pip for Python 3.11 using get-pip.py
RUN wget https://bootstrap.pypa.io/get-pip.py \
    && python3.11 get-pip.py \
    && rm get-pip.py

# Set Python 3.11 as the default python and pip
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.11 1 \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1 \
    && ln -sf /usr/local/bin/pip3.11 /usr/bin/pip

# Upgrade pip to latest version
RUN python -m pip install --upgrade pip

# Set working directory
WORKDIR /workspace

# Clone the Dhaka repository and install
RUN git clone https://github.com/BoevaLab/Dhaka-implementation.git \
    && cd Dhaka-implementation \
    && pip install --no-cache-dir . scanpy

# Set the working directory back to workspace
WORKDIR /workspace
