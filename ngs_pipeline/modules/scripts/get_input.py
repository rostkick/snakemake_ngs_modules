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

    # Get calling mode from config
    calling_mode = config.get('calling_mode', 'both')  # 'individual', 'joint', 'both'

    if ngs.GRM:
        
        archive_bams = [
            f"archive/{config['run']}/bam/{sample}.final.bam"
            for sample in ngs.SAMPLES
        ]
        archive_bai = [
            f"archive/{config['run']}/bam/{sample}.final.bam.bai"
            for sample in ngs.SAMPLES
        ]
        germline_inputs += archive_bams

        # Joint calling (cohort-level analysis)
        if calling_mode in ['joint', 'both']:
            germline_inputs_cohort = [
                f"results/{config['run']}/germline/vcf/cohort.annotated.vcf.gz"
            ]
            germline_inputs += germline_inputs_cohort


        # Individual calling (per-sample analysis)
        if calling_mode in ['individual', 'both']:
            germline_inputs_individual = [
                f"results/{config['run']}/germline/xlsx/individual.{config['run']}.germline.results.xlsx"
            ]
            germline_inputs += germline_inputs_individual


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

    # Quality metrics for WES samples
    if config['ngs_type'] == 'WES':
        metrics = [
            f"results/{config['run']}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]

    input_files = germline_inputs + somatic_inputs + metrics
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
    
    # Get calling mode from config
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
    Get quality metrics targets for WES samples.
    
    Args:
        ngs: NGS data object containing sample information
        
    Returns:
        list: List of metrics target files
    """
    if config['ngs_type'] == 'WES':
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
    # Temporarily override config for joint-only mode
    original_mode = config.get('calling_mode', 'both')
    config['calling_mode'] = 'joint'
    
    result = get_final_inputs(ngs)
    
    # Restore original config
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
    # Temporarily override config for individual-only mode
    original_mode = config.get('calling_mode', 'both')
    config['calling_mode'] = 'individual'
    
    result = get_final_inputs(ngs)
    
    # Restore original config
    config['calling_mode'] = original_mode
    
    return result