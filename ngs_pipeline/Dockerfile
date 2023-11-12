## ===================
FROM ngs.deps AS ngs.final

WORKDIR /usr/ngs

COPY ngs/snakefile.smk ./pipeline/snakefile.smk
COPY ngs/configure_run.yml ./pipeline/configure_run.yml
COPY ngs/configure_tools.yml ./pipeline/configure_tools.yml
COPY ngs/configure_reference_37.yml ./pipeline/configure_reference_37.yml
COPY ngs/configure_reference_38.yml ./pipeline/configure_reference_38.yml
COPY ngs/modules ./pipeline/modules/
