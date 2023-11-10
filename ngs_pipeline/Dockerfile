## ===================
FROM ngs.deps AS ngs.final

WORKDIR /usr/ngs

COPY ngs/snakefile.smk ./pipeline/snakefile.smk
COPY ngs/configure.yml ./pipeline/configure.yml
COPY ngs/modules ./pipeline/modules/
