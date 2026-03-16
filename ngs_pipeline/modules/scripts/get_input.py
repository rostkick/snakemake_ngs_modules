from itertools import product
from snakemake.workflow import config


def get_final_inputs(ngs):
    germline_inputs = []
    somatic_inputs = []
    metrics = []

    calling_mode = config.get('calling_mode', 'joint')
    ngs_type = config['ngs_type']
    run = config['run']

    _panel_names = [
        p['name'] for p in config['references'].get('wes_gene_panels', [])
    ] if ngs_type == 'WES' else []

    if ngs.GRM:

        # Archive dedup BAM (DeepVariant input, no BQSR) and BQSR BAM (Mutect2 input)
        germline_inputs += [
            f"archive/{run}/bam/{sample}.dedup.bam"
            for sample in ngs.SAMPLES
        ]
        germline_inputs += [
            f"archive/{run}/bam/{sample}.final.bqsr.bam"
            for sample in ngs.SAMPLES
        ]

        # Archive per-sample DeepVariant VCF + gVCF (r0_2, individual mode only)
        if calling_mode == 'individual':
            germline_inputs += [
                f"archive/{run}/germline/vcf/{sample}.vcf.gz"
                for sample in ngs.GRM_SAMPLES
            ]

        # Joint calling — cohort-level annotated VCF (r8_1)
        if calling_mode == 'joint':
            germline_inputs += [
                f"results/{run}/germline/vcf/cohort.annotated.vcf.gz"
            ]

        if ngs_type == 'WES':
            # Per-panel xlsx (r9_8)
            germline_inputs += [
                f"results/{run}/germline/xlsx/panel.{panel_name}.{run}.germline.results.xlsx"
                for panel_name in _panel_names
            ]
            # Per-sample full exome xlsx (r9_9)
            germline_inputs += [
                f"results/{run}/germline/xlsx/wes_clinical.{run}.germline.results.xlsx"
            ]
            # Archive panel TSVs (r0_4)
            germline_inputs += [
                f"archive/{run}/germline/tsv/{sample}.{panel_name}.tsv"
                for sample in ngs.GRM_SAMPLES
                for panel_name in _panel_names
            ]
            # Archive full WES TSVs (r0_5)
            germline_inputs += [
                f"archive/{run}/germline/tsv/{sample}.wes_clinical.tsv"
                for sample in ngs.GRM_SAMPLES
            ]
            # Archive raw unfiltered TSVs (r0_9)
            germline_inputs += [
                f"archive/{run}/germline/tsv/{sample}.raw.tsv"
                for sample in ngs.GRM_SAMPLES
            ]
            # Archive panel xlsx (r0_10)
            germline_inputs += [
                f"archive/{run}/germline/xlsx/panel.{panel_name}.{run}.germline.results.xlsx"
                for panel_name in _panel_names
            ]
            # Archive wes_clinical xlsx (r0_11)
            germline_inputs += [
                f"archive/{run}/germline/xlsx/wes_clinical.{run}.germline.results.xlsx"
            ]
            # Archive gVCFs for future joint calling (r0_12)
            germline_inputs += [
                f"archive/{run}/germline/gvcf/{sample}.gvcf.gz"
                for sample in ngs.GRM_SAMPLES
            ]
            # Archive annotated cohort VCF (r0_13)
            germline_inputs += [
                f"archive/{run}/germline/vcf/cohort.annotated.vcf.gz"
            ]

        else:
            # panel ngs_type and WGS — single cohort xlsx (r9_10) + archived TSVs (r0_3)
            germline_inputs += [
                f"results/{run}/germline/xlsx/individual.{run}.germline.results.xlsx",
                f"archive/{run}/germline/xlsx/individual.{run}.germline.results.xlsx",
            ]
            germline_inputs += [
                f"archive/{run}/germline/tsv/{sample}.tsv"
                for sample in ngs.GRM_SAMPLES
            ]

        # Somatic calling — paired germline vs tumour
        if ngs.TMR:
            somatic_inputs += [
                f"results/{run}/somatic/{patient}/somatic_annotated.vcf.gz"
                for patient in ngs.GRM_VS_TMR_PATIENTS
            ]

    # Tumour-only samples
    if len(ngs.ONLY_TMR_PATIENTS) > 0:
        somatic_inputs += [
            f"results/{run}/somatic/{patient}/somatic_annotated_tonly.vcf.gz"
            for patient in ngs.ONLY_TMR_PATIENTS
        ]

    # Hybrid-capture quality metrics for WES and panel samples
    if ngs_type in ['WES', 'panel']:
        metrics += [
            f"results/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]
        metrics += [
            f"archive/{run}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]

    # FASTQ checksums (r0_7)
    provenance = [
        f"results/{run}/provenance/fastq_checksums.md5"
    ]

    return germline_inputs + somatic_inputs + metrics + provenance


def get_germline_calling_targets(ngs):
    calling_mode = config.get('calling_mode', 'joint')
    ngs_type = config['ngs_type']
    run = config['run']
    targets = []

    _panel_names = [
        p['name'] for p in config['references'].get('wes_gene_panels', [])
    ] if ngs_type == 'WES' else []

    if ngs.GRM:
        if calling_mode == 'joint':
            targets.append(f"results/{run}/germline/vcf/cohort.annotated.vcf.gz")
        if ngs_type == 'WES':
            targets += [
                f"results/{run}/germline/xlsx/panel.{panel_name}.{run}.germline.results.xlsx"
                for panel_name in _panel_names
            ]
            targets += [
                f"results/{run}/germline/xlsx/wes_clinical.{run}.germline.results.xlsx"
            ]
        else:
            targets.append(
                f"results/{run}/germline/xlsx/individual.{run}.germline.results.xlsx"
            )

    return targets


def get_somatic_calling_targets(ngs):
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
    if config['ngs_type'] in ['WES', 'panel']:
        return [
            f"results/{config['run']}/bam/hs_metrics/{sample}.hs_metrics.tsv"
            for sample in ngs.SAMPLES
        ]
    return []