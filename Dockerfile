# # First stage: Build the initial image with the all heavy dependencies
# FROM continuumio/miniconda3 AS deps

# # Copy the environment.yml file & tools directory into the container
# COPY deps/environment.yml /ngs_pipeline/environment.yml
# COPY deps/tools /ngs_pipeline/tools

# RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
# 					conda init bash && \
# 					conda activate base && \
# 					conda install --yes conda-forge::libarchive && \
# 					conda install --yes conda-forge::mamba && \
# 					mamba env create -f /ngs_pipeline/environment.yml"
# # RUN conda init bash
# # RUN conda activate base

# # # Get mamba deps
# # RUN conda install conda-forge::libarchive

# # # Get mamba for further installation acceleration 
# # RUN conda install --yes conda-forge::mamba
# # RUN mamba env create -f /ngs_pipeline/environment.yml

# # # Init Conda & clean up cache to reduce image size
# # RUN /bin/bash -c "source /root/.bashrc && mamba activate base && mamba env create -f /ngs_pipeline/environment.yml" && \
# #     mamba env update -n base --file /ngs_pipeline/environment.yml && \
# #     mamba clean -a -y

# # ============================================================
# # Second stage: Build the final image with the code
# FROM deps AS code

# # Copy the project code into the container
# COPY /ngs_pipeline /ngs_pipeline

# # Set the working directory
# WORKDIR /ngs_pipeline

# CMD ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate smk && snakemake --snakefile Snakefile --cores 4; tail -f /dev/null"]
# # SHELL [ "conda", "activate", "smk"]


# =====================
# Multi-stage build

# Stage 1: Dependencies
FROM continuumio/miniconda3 as deps

# Установка зависимостей в контейнере
WORKDIR /ngs_pipeline

# Deps for singularity
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

# Копируем файл с зависимостями и инструменты
COPY deps/environment.yml .
COPY deps/tools ./tools

# Mamba
RUN conda install conda-forge::libarchive
RUN conda install --yes conda-forge::mamba

# Создаем окружение из environment.yml
RUN mamba env create -f environment.yml

# Подготавливаем окружение, чтобы оно использовалось по умолчанию
SHELL ["conda", "run", "-n", "smk", "/bin/bash", "-c"]
RUN conda init bash

# Stage 2: Application
# FROM continuumio/miniconda3 as code
FROM deps as code

# Копируем только окружение из первого этапа
# COPY --from=deps /opt/conda /opt/conda

# Устанавливаем рабочую директорию
WORKDIR /ngs_pipeline

# Копируем весь код проекта
COPY ngs_pipeline app

# Активируем окружение по умолчанию
SHELL ["conda", "run", "-n", "smk", "/bin/bash", "-c"]

# Устанавливаем необходимые настройки и активируем окружение
RUN echo "conda activate smk" >> ~/.bashrc

# Команда по умолчанию
CMD ["bash"]
