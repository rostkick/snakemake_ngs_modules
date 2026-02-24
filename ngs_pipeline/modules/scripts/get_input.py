from itertools import product
from snakemake.workflow import config


def get_final_inputs(ngs):
    """
    Return a list of input files for all samples based on calling strategy.

    Args:
        ngs: NGS data object containing sample information

    Returns:
        list: List of target files based on configuration
    """
    germline_inputs = []
    somatic_inputs = []
    metrics = []

    calling_mode = config.get('calling_mode', 'both')

    if ngs.GRM:

        # Archive BAM files (rule r0_1_archive_bams)
        archive_bams = [
            f"archive/{config['run']}/bam/{sample}.final.bam"
            for sample in ngs.SAMPLES
        ]
        germline_inputs += archive_bams

        # Archive per-sample DeepVariant VCF + gVCF (rule r0_2_archive_vcfs)
        if calling_mode in ['individual', 'both']:
            archive_vcfs = [
                f"archive/{config['run']}/germline/vcf/{sample}.vcf.gz"
                for sample in ngs.GRM_SAMPLES
            ]
            germline_inputs += archive_vcfs

        # Joint calling — cohort-level annotated VCF
        if calling_mode in ['joint', 'both']:
            germline_inputs += [
                f"results/{config['run']}/germline/vcf/cohort.annotated.vcf.gz"
            ]

        # Individual calling — per-sample XLSX + archived TSV
        if calling_mode in ['individual', 'both']:
            germline_inputs += [
                f"results/{config['run']}/germline/xlsx/individual.{config['run']}.germline.results.xlsx"
            ]

            # Archive XLSX (rule r0_4_archive_xlsx)
            germline_inputs += [
                f"archive/{config['run']}/germline/xlsx/individual.{config['run']}.germline.results.xlsx"
            ]

            # Archive TSV files (rule r0_3_archive_tsv_individual)
            germline_inputs += [
                f"archive/{config['run']}/germline/tsv/{sample}.tsv"
                for sample in ngs.GRM_SAMPLES
            ]

        # Somatic calling — paired germline vs tumour
        if ngs.TMR:
            somatic_inputs += [
                f"results/{config['run']}/somatic/{patient}/somatic_annotated.vcf.gz"
                for patient in ngs.GRM_VS_TMR_PATIENTS
            ]

    # Tumour-only samples
    if len(ngs.ONLY_TMR_PATIENTS) > 0:
        somatic_inputs += [
            f"results/{config['run']}/somatic/{patient}/somatic_annotated_tonly.vcf.gz"
            for patient in ngs.ONLY_TMR_PATIENTS
        ]

    # Hybrid-capture quality metrics for WES and panel samples
    if config['ngs_type'] in ['WES', 'panel']:
        metrics = [
            f"results/{config['run']}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]
        # Archive HS metrics (rule r0_7_archive_hs_metrics)
        metrics += [
            f"archive/{config['run']}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]

    # FASTQ checksums — always produced for full input traceability
    provenance = [
        f"results/{config['run']}/provenance/fastq_checksums.md5"
    ]

    return germline_inputs + somatic_inputs + metrics + provenance


def get_germline_calling_targets(ngs):
    """
    Get only germline calling targets based on configuration.

    Args:
        ngs: NGS data object containing sample information

    Returns:
        list: List of germline target files
    """
    targets = []

    calling_mode = config.get('calling_mode', 'both')

    if ngs.GRM:
        if calling_mode in ['joint', 'both']:
            targets.append(f"results/{config['run']}/germline/vcf/cohort.annotated.vcf.gz")

        if calling_mode in ['individual', 'both']:
            targets.append(
                f"results/{config['run']}/germline/xlsx/individual.{config['run']}.germline.results.xlsx"
            )

    return targets


def get_somatic_calling_targets(ngs):
    """
    Get only somatic calling targets.

    Args:
        ngs: NGS data object containing sample information

    Returns:
        list: List of somatic target files
    """
    targets = []

    if ngs.GRM and ngs.TMR:
        targets += [
            f"results/{config['run']}/somatic/{patient}/somatic_annotated.vcf.gz"
            for patient in ngs.GRM_VS_TMR_PATIENTS
        ]

    if len(ngs.ONLY_TMR_PATIENTS) > 0:
        targets += [
            f"results/{config['run']}/somatic/{patient}/somatic_annotated_tonly.vcf.gz"
            for patient in ngs.ONLY_TMR_PATIENTS
        ]

    return targets


def get_metrics_targets(ngs):
    """
    Get quality metrics targets for WES and panel samples.

    Args:
        ngs: NGS data object containing sample information

    Returns:
        list: List of metrics target files
    """
    if config['ngs_type'] in ['WES', 'panel']:
        return [
            f"results/{config['run']}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]
    return []


def get_joint_only_inputs(ngs):
    """
    Get inputs for joint calling only mode.

    Args:
        ngs: NGS data object containing sample information

    Returns:
        list: List of target files for joint calling only
    """
    original_mode = config.get('calling_mode', 'both')
    config['calling_mode'] = 'joint'

    result = get_final_inputs(ngs)

    config['calling_mode'] = original_mode

    return result


def get_individual_only_inputs(ngs):
    """
    Get inputs for individual calling only mode.

    Args:
        ngs: NGS data object containing sample information

    Returns:
        list: List of target files for individual calling only
    """
    original_mode = config.get('calling_mode', 'both')
    config['calling_mode'] = 'individual'

    result = get_final_inputs(ngs)

    config['calling_mode'] = original_mode

    return result