# =====================
# Multi-stage build

# Stage 1: Dependencies
FROM continuumio/miniconda3 as deps

# Set working directory
WORKDIR /ngs_pipeline

# Install system dependencies including Java
RUN apt-get update && apt-get install -y \
		build-essential \
		libssl-dev \
		uuid-dev \
		libgpgme11-dev \
		squashfs-tools \
		libseccomp-dev \
		gperf \
		wget \
		pkg-config \
		git \
		rsync \
		default-jdk \
		locales \
		ca-certificates && \
	rm -rf /var/lib/apt/lists/*

# Configure locale
RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen en_US.UTF-8 && \
    update-locale LANG=en_US.UTF-8

ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Copy environment file and tools
COPY deps/environment.yml .
COPY deps/tools ./tools

# go installation
RUN cd /usr/local && \
		wget https://golang.org/dl/go1.16.5.linux-amd64.tar.gz && \
		rm -rf /usr/local/go && tar -C /usr/local -xzf go1.16.5.linux-amd64.tar.gz

# libseccomp installation - needed to compile singularity
RUN wget https://github.com/seccomp/libseccomp/releases/download/v2.5.1/libseccomp-2.5.1.tar.gz && \
tar zxvf libseccomp-2.5.1.tar.gz && rm libseccomp-2.5.1.tar.gz && \
cd libseccomp-2.5.1 && \
./configure && make && make install && cd ~
RUN apt update && apt-get upgrade -y && apt-get install libseccomp-dev -y
ENV PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig/:$PKG_CONFIG_PATH

# singularity 
RUN export VERSION=3.7.4 && \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-${VERSION}.tar.gz && \
    tar -xzf singularity-${VERSION}.tar.gz && rm singularity-${VERSION}.tar.gz

RUN cd singularity && \
    export PATH=$PATH:/usr/local/go/bin && \
    /bin/bash -c "source ~/.profile" && \
    ./mconfig && \ 
    make BUILDTAGS="" -C ./builddir && \
    make -C ./builddir install

# Install Mamba
RUN conda config --set solver libmamba && \
    conda install --yes conda-forge::mamba

# Create conda environment from environment.yml
RUN mamba env create --yes -f environment.yml
RUN mamba install -c conda-forge gsl -y

# Remove conda's openjdk and reinstall picard without openjdk dependency
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate smk && \
    conda remove --force -y openjdk || true && \
    mamba install -c bioconda --yes picard --no-deps && \
    conda clean -a -y"

# RUN mamba install -c conda-forge -c bioconda bcftools
# RUN ln -sf /opt/conda/lib/libgsl.so.25 /opt/conda/envs/smk/lib/libgsl.so.25 && ln -sf /opt/conda/lib/libgslcblas.so.0 /opt/conda/envs/smk/lib/libgslcblas.so.0

# Initialize conda for bash
RUN conda init bash

# Stage 2: Application
FROM deps as code

# Set working directory
WORKDIR /ngs_pipeline

# Copy project code
COPY ngs_pipeline app

# Find and set JAVA_HOME dynamically
RUN JAVA_PATH=$(update-alternatives --query java | grep 'Value:' | cut -d' ' -f2) && \
    JAVA_HOME=$(dirname $(dirname $JAVA_PATH)) && \
    echo "export JAVA_HOME=$JAVA_HOME" >> ~/.bashrc && \
    echo "export PATH=\$JAVA_HOME/bin:/opt/conda/envs/smk/bin:\$PATH" >> ~/.bashrc

ENV PATH=/opt/conda/envs/smk/bin:$PATH

# Activate conda environment by default
RUN echo "conda activate smk" >> ~/.bashrc

# Default command
CMD ["bash"]