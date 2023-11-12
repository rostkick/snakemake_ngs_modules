# Use the official image as a parent image.
FROM ubuntu:22.04 AS ngs.base

# Set the working directory.
USER root
WORKDIR /usr/ngs
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y wget tar htop vim

# Setting up conda and apt-based core libs
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /usr/local/conda
ENV PATH=/usr/local/conda/bin:$PATH
RUN apt-get install -y apt-transport-https 
RUN apt-get update -y && apt-get install -y software-properties-common && apt-get install -y build-essential && apt-get install -y libz-dev
RUN conda update -n base -c defaults conda -y
RUN conda install mamba -c conda-forge -y

# Final setup of environment
RUN mamba install -c anaconda python>=3.10
RUN mamba install -c bioconda -c conda-forge snakemake -y

FROM ngs.base AS ngs.deps

# Some NGS software - species detection and right-wing
RUN mamba install -c bioconda bwa-mem2 samtools gatk4 seqkit bedtools pandas -y

## VEP install
# deps for singularity

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
    git

# go installation
RUN cd /usr/local && \
    wget https://golang.org/dl/go1.16.5.linux-amd64.tar.gz && \
    rm -rf /usr/local/go && tar -C /usr/local -xzf go1.16.5.linux-amd64.tar.gz

RUN cd /usr/ngs

# libseccomp installation - needed to compile singularity
RUN wget https://github.com/seccomp/libseccomp/releases/download/v2.5.1/libseccomp-2.5.1.tar.gz && \
    tar zxvf libseccomp-2.5.1.tar.gz && rm libseccomp-2.5.1.tar.gz && \
    cd libseccomp-2.5.1 && \
    ./configure && make && make install && cd /usr/ngs

#RUN echo 'deb http://deb.debian.org/debian buster-backports main' >> /etc/apt/sources.list
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

COPY tools /usr/ngs/tools
