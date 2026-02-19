"""Compute MD5 checksums for all input FASTQ files in a run.
Ensures full traceability of input data.

Snakemake script interface:
    snakemake.output.md5   - output checksum file
    snakemake.config       - pipeline config (grm_dir, tmr_dir)
    snakemake.log[0]       - log file
"""
import os
import hashlib


def compute_fastq_checksums(config, output_md5, log_path):
    fastq_dirs = []
    grm_dir = config.get('grm_dir', '')
    tmr_dir = config.get('tmr_dir', '')
    if grm_dir:
        fastq_dirs.append(grm_dir)
    if tmr_dir:
        fastq_dirs.append(tmr_dir)

    checksums = []
    for fq_dir in fastq_dirs:
        if not os.path.isdir(fq_dir):
            continue
        for fname in sorted(os.listdir(fq_dir)):
            if not any(fname.endswith(ext) for ext in ('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
                continue
            fpath = os.path.join(fq_dir, fname)
            # Follow symlinks to compute checksum of actual file
            real_path = os.path.realpath(fpath)
            md5 = hashlib.md5()
            with open(real_path, 'rb') as f:
                for chunk in iter(lambda: f.read(8192), b''):
                    md5.update(chunk)
            checksums.append(f"{md5.hexdigest()}  {real_path}  # link: {fpath}")

    os.makedirs(os.path.dirname(output_md5), exist_ok=True)
    with open(output_md5, 'w') as out:
        for line in checksums:
            out.write(line + '\n')
    with open(log_path, 'w') as logf:
        logf.write(f"Computed {len(checksums)} FASTQ checksums\n")


compute_fastq_checksums(snakemake.config, snakemake.output.md5, str(snakemake.log))
