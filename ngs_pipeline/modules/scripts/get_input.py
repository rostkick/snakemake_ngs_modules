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
        
        # Archive VCF files (rule r0_2_archive_vcfs)
        if calling_mode in ['individual', 'both']:
            archive_vcfs = [
                f"archive/{config['run']}/germline/vcf/{sample}.vcf.gz"
                for sample in ngs.GRM_SAMPLES
            ]
            germline_inputs += archive_vcfs

        # Joint calling (cohort-level analysis)
        if calling_mode in ['joint', 'both']:
            germline_inputs_cohort = [
                f"results/{config['run']}/germline/vcf/cohort.annotated.vcf.gz"
            ]
            germline_inputs += germline_inputs_cohort

        # Archive GATK HC GVCF (always when joint/both — valuable for future joint calling)
        if calling_mode in ['joint', 'both']:
            gvcf_suffix = '.wgs.gatk.gvcf.gz' if config['ngs_type'] == 'WGS' else '.gatk.gvcf.gz'
            archive_gatk_gvcfs = [
                f"archive/{config['run']}/germline/gvcf/{sample}{gvcf_suffix}"
                for sample in ngs.GRM_SAMPLES
            ]
            germline_inputs += archive_gatk_gvcfs

        # Individual calling (per-sample analysis)
        if calling_mode in ['individual', 'both']:
            # Final XLSX (remains in results)
            germline_inputs_individual = [
                f"results/{config['run']}/germline/xlsx/individual.{config['run']}.germline.results.xlsx"
            ]
            germline_inputs += germline_inputs_individual
            
            # Archive XLSX (rule r0_4_archive_xlsx)
            archive_xlsx = [
                f"archive/{config['run']}/germline/xlsx/individual.{config['run']}.germline.results.xlsx"
            ]
            germline_inputs += archive_xlsx
            
            # Archive TSV files (rule r0_3_archive_tsv_individual)
            archive_tsvs = [
                f"archive/{config['run']}/germline/tsv/{sample}.tsv"
                for sample in ngs.GRM_SAMPLES
            ]
            germline_inputs += archive_tsvs

        # Somatic calling (if tumor samples are present)
        if ngs.TMR:
            somatic_inputs += [
                f"results/{config['run']}/somatic/{patient}/somatic_annotated.vcf.gz"
                for patient in ngs.GRM_VS_TMR_PATIENTS
            ]

    # Tumor-only samples processing
    if len(ngs.ONLY_TMR_PATIENTS) > 0:
        somatic_inputs += [
            f"results/{config['run']}/somatic/{patient}/somatic_annotated_tonly.vcf.gz"
            for patient in ngs.ONLY_TMR_PATIENTS
        ]

    # Quality metrics for WES and panel samples
    if config['ngs_type'] in ['WES', 'panel']:
        metrics = [
            f"results/{config['run']}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]
        # Archive HS metrics
        archive_hs = [
            f"archive/{config['run']}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]
        metrics += archive_hs

    # FASTQ checksums (always — for full traceability)
    provenance = [
        f"results/{config['run']}/provenance/fastq_checksums.md5"
    ]

    input_files = germline_inputs + somatic_inputs + metrics + provenance
    return input_files


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
        # Joint calling targets
        if calling_mode in ['joint', 'both']:
            targets.append(f"results/{config['run']}/germline/vcf/cohort.annotated.vcf.gz")
        
        # Individual calling targets
        if calling_mode in ['individual', 'both']:
            targets.append(f"results/{config['run']}/germline/xlsx/individual.{config['run']}.germline.results.xlsx")
    
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
    
    # Paired germline vs tumor samples
    if ngs.GRM and ngs.TMR:
        targets += [
            f"results/{config['run']}/somatic/{patient}/somatic_annotated.vcf.gz"
            for patient in ngs.GRM_VS_TMR_PATIENTS
        ]
    
    # Tumor-only samples
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