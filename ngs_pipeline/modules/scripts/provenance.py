"""
provenance.py — Collect reproducibility metadata for each pipeline run.

Called from snakefile onstart: block.
Writes all provenance data into results/{run}/provenance/.
"""

import os
import json
import subprocess
import platform
from datetime import datetime
from snakemake.workflow import config


def _run_cmd(cmd: str, default: str = "N/A") -> str:
    """Run a shell command and return stripped stdout, or default on failure."""
    try:
        return subprocess.check_output(
            cmd, shell=True, stderr=subprocess.DEVNULL, text=True
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return default


def check_git_clean(force: bool = False) -> dict:
    """
    Check if the pipeline repository is clean (all changes committed).

    Raises SystemExit if dirty and force is False.
    Returns dict with git status info.
    """
    git_sha = _run_cmd("git rev-parse HEAD")
    git_branch = _run_cmd("git rev-parse --abbrev-ref HEAD")
    git_dirty = _run_cmd("git status --porcelain")

    is_dirty = len(git_dirty) > 0

    info = {
        "git_sha": git_sha,
        "git_branch": git_branch,
        "git_dirty": is_dirty,
        "git_dirty_files": git_dirty if is_dirty else "",
    }

    if is_dirty and not force:
        print("\n" + "=" * 70)
        print("  ERROR: Pipeline repository has uncommitted changes!")
        print("=" * 70)
        print(f"  Branch : {git_branch}")
        print(f"  SHA    : {git_sha}")
        print(f"\n  Dirty files:")
        for line in git_dirty.splitlines():
            print(f"    {line}")
        print(f"\n  Options:")
        print(f"    1. Commit changes:   git add -A && git commit -m 'pre-run snapshot'")
        print(f"    2. Force (not recommended): snakemake ... --config force_dirty=True")
        print("=" * 70 + "\n")
        raise SystemExit(
            "Aborting: uncommitted changes. Use --force-dirty to override."
        )

    if is_dirty and force:
        print(f"\n  WARNING: Running with dirty repository (--force-dirty).")
        print(f"           Results may not be fully reproducible.\n")

    return info


def collect_tool_versions() -> dict:
    """Collect versions of all bioinformatics tools used by the pipeline."""
    tools = config.get("tools", {})
    versions = {}

    # samtools
    samtools = tools.get("samtools", "samtools")
    versions["samtools"] = _run_cmd(f"{samtools} --version | head -1")

    # bwa-mem2
    bwa = tools.get("bwa_mem2", "bwa-mem2")
    versions["bwa-mem2"] = _run_cmd(f"{bwa} version 2>&1 | head -1")

    # gatk
    gatk = tools.get("gatk", "gatk")
    versions["gatk"] = _run_cmd(f"{gatk} --version 2>&1 | grep 'GATK' | head -1")

    # picard
    picard = tools.get("picard", "picard")
    versions["picard"] = _run_cmd(f"{picard} MarkDuplicates --version 2>&1 | head -1")

    # bcftools
    bcftools = tools.get("bcftools", "bcftools")
    versions["bcftools"] = _run_cmd(f"{bcftools} --version | head -1")

    # snakemake
    versions["snakemake"] = _run_cmd("snakemake --version")

    # python
    versions["python"] = platform.python_version()

    # singularity (for DeepVariant/VEP containers)
    singularity = tools.get("singularity", "singularity")
    versions["singularity"] = _run_cmd(f"{singularity} --version 2>&1")

    # Container image digests (SIF files)
    dv_sif = tools.get("deepvariant", "")
    if dv_sif and os.path.exists(dv_sif):
        versions["deepvariant_sif"] = dv_sif
        versions["deepvariant_sif_md5"] = _run_cmd(f"md5sum {dv_sif} | cut -d' ' -f1")
        versions["deepvariant_sif_size"] = os.path.getsize(dv_sif)

    vep_conf = tools.get("vep", {})
    vep_sif = vep_conf.get("path", "") if isinstance(vep_conf, dict) else ""
    if vep_sif and os.path.exists(vep_sif):
        versions["vep_sif"] = vep_sif
        versions["vep_sif_md5"] = _run_cmd(f"md5sum {vep_sif} | cut -d' ' -f1")
        versions["vep_sif_size"] = os.path.getsize(vep_sif)

    return versions


def collect_environment_info() -> dict:
    """Collect execution environment details."""
    return {
        "hostname": platform.node(),
        "os": f"{platform.system()} {platform.release()}",
        "architecture": platform.machine(),
        "cpu_count": os.cpu_count(),
        "date_utc": datetime.utcnow().isoformat() + "Z",
        "date_local": datetime.now().isoformat(),
        "user": os.environ.get("USER", "unknown"),
        "cwd": os.getcwd(),
    }


def collect_reference_info() -> dict:
    """Collect reference file paths and sizes (md5 of genome is too slow)."""
    refs = config.get("references", {})
    info = {}
    for key, path in refs.items():
        if isinstance(path, str) and os.path.exists(path):
            stat = os.stat(path)
            info[key] = {
                "path": path,
                "size_bytes": stat.st_size,
                "mtime": datetime.fromtimestamp(stat.st_mtime).isoformat(),
            }
        elif isinstance(path, dict):
            # nested config like vep_plugins_data
            info[key] = path
        else:
            info[key] = {"path": str(path), "exists": False}
    return info


def collect_conda_env() -> str:
    """Collect conda environment package list."""
    return _run_cmd("conda list 2>/dev/null", default="conda not available")


def dump_merged_config() -> dict:
    """Return the fully merged Snakemake config (all configfiles merged)."""
    # Filter out non-serializable items
    serializable = {}
    for k, v in config.items():
        try:
            json.dumps(v)
            serializable[k] = v
        except (TypeError, ValueError):
            serializable[k] = str(v)
    return serializable


def write_provenance(run_name: str, force_dirty: bool = False):
    """
    Main entry point: collect all provenance data and write to disk.

    Called from snakefile onstart: block.
    """
    prov_dir = f"results/{run_name}/provenance"
    os.makedirs(prov_dir, exist_ok=True)

    # 1. Git check (strict — will exit if dirty and not forced)
    git_info = check_git_clean(force=force_dirty)

    # 2. Tool versions
    tool_versions = collect_tool_versions()

    # 3. Execution environment
    env_info = collect_environment_info()

    # 4. Merged config (exact parameters of this run)
    merged_config = dump_merged_config()

    # 5. Reference files info
    ref_info = collect_reference_info()

    # Assemble provenance record
    provenance = {
        "pipeline_version": {
            "git_sha": git_info["git_sha"],
            "git_branch": git_info["git_branch"],
            "git_dirty": git_info["git_dirty"],
            "git_dirty_files": git_info["git_dirty_files"],
            "force_dirty": force_dirty,
        },
        "tool_versions": tool_versions,
        "environment": env_info,
        "references": ref_info,
        "run_name": run_name,
    }

    # Write provenance JSON
    prov_file = os.path.join(prov_dir, "provenance.json")
    with open(prov_file, "w") as f:
        json.dump(provenance, f, indent=2, default=str)
    print(f"  Provenance written: {prov_file}")

    # Write merged config YAML
    config_file = os.path.join(prov_dir, "merged_config.json")
    with open(config_file, "w") as f:
        json.dump(merged_config, f, indent=2, default=str)
    print(f"  Merged config written: {config_file}")

    # Write conda environment snapshot
    conda_snapshot = collect_conda_env()
    conda_file = os.path.join(prov_dir, "conda_packages.txt")
    with open(conda_file, "w") as f:
        f.write(conda_snapshot)
    print(f"  Conda snapshot written: {conda_file}")

    # Write DAG description (rule list)
    dag_file = os.path.join(prov_dir, "git_sha.txt")
    with open(dag_file, "w") as f:
        f.write(git_info["git_sha"] + "\n")
    print(f"  Git SHA written: {dag_file}")

    print(f"\n  Provenance collection complete for run '{run_name}'")
    print(f"  SHA: {git_info['git_sha']}")
    print(f"  Dirty: {git_info['git_dirty']}\n")
