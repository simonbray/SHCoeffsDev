# Debian base image
FROM debian:bullseye-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH

# Install base dependencies and micromamba
RUN apt-get update && apt-get install -y \
    curl \
    wget \
    bzip2 \
    git \
    build-essential \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

# Install Micromamba
RUN curl -L https://micromamba.snakepit.net/api/micromamba/linux-64/latest | tar -xvj -C /usr/local/bin --strip-components=1 bin/micromamba

# Create the environment and install packages
RUN micromamba create -y -p $MAMBA_ROOT_PREFIX/envs/pyenv \
    python=2.7.15 \
    rdkit=2018.09.3 \
    keras=2.2.4 \
    tensorflow=1.10.0 \
    numpy=1.16.2 \
    xarray=0.11.3 \
    ##--channel=rdkit \
    --channel=conda-forge \
    --channel=defaults

RUN micromamba create -y -p $MAMBA_ROOT_PREFIX/envs/lirasearch \
    python=3.11 \
    tqdm=4.67.1 \
    pyignite=0.6.1 \
    scipy=1.15.2 \
    colorama=0.4.6 \
    numpy=2.2.6 \
    --channel=conda-forge \
    --channel=defaults

# Activate environment by default
ENV PATH=$MAMBA_ROOT_PREFIX/envs/pyenv/bin:$PATH

# Install Node.js and npm v18
RUN curl -fsSL https://deb.nodesource.com/setup_18.x | bash - && \
    apt-get update && \
    apt-get install -y nodejs

# Install Julia (in /opt/, not as root)
ENV JULIA_VERSION=1.11.4
ENV JULIA_DEPOT_PATH=/opt/.julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.11/julia-$JULIA_VERSION-linux-x86_64.tar.gz && \
    tar -xvzf julia-$JULIA_VERSION-linux-x86_64.tar.gz && \
    mv julia-$JULIA_VERSION /opt/julia && \
    ln -s /opt/julia/bin/julia /usr/local/bin/julia && \
    rm julia-$JULIA_VERSION-linux-x86_64.tar.gz

# Install Julia packages
ENV JULIA_CPU_TARGET=x86_64;haswell;skylake;skylake-avx512;tigerlake
RUN julia -e 'using Pkg; Pkg.add.(["Glob", "Printf", "ArgParse", "TimerOutputs", "GeometryBasics", "ImplicitBVH", "Distributions"])'
RUN julia -e 'using Pkg; Pkg.precompile()'
#RUN cp -r /root/.julia/environments/v1.11 /opt/juliaenv/

# Create and switch to GitRepos directory
RUN mkdir /opt/GitRepos
WORKDIR /opt/GitRepos

# Clone git repos
RUN git clone https://github.com/AstexUK/ESP_DNN.git
RUN git clone https://github.com/AstexUK/esp-surface-generator.git
#ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN git clone https://github.com/simonbray/SHCoeffsDev.git

# Install from the git repos
RUN cd ESP_DNN && python setup.py install && cd ..
RUN cd esp-surface-generator/ && npm install && cd .. && ln -s /opt/GitRepos/esp-surface-generator/cli.js /usr/bin/esp-surface-generator
RUN ln -s /opt/GitRepos/SHCoeffsDev/sh_coeffs.jl /usr/bin/sh-coeff-calculator && ln -s /opt/GitRepos/SHCoeffsDev/Geotools.jl /usr/bin/ 
RUN chmod +x /opt/GitRepos/SHCoeffsDev/sh_coeffs.jl
WORKDIR /

# Barely achieves anything, but clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN micromamba clean --all --yes

# Set default shell to bash
SHELL ["/bin/bash", "-c"]

# Default command
CMD ["bash"]

