import fcntl
import os
import shutil
import time
from pathlib import Path


def _as_single(value, name):
    if isinstance(value, (list, tuple)):
        if len(value) != 1:
            raise ValueError(f"Expected a single path for {name}, got {len(value)}")
        return str(value[0])
    return str(value)


def _log(log_path, message):
    with open(log_path, "a", encoding="utf-8") as handle:
        handle.write(message + "\n")


def _write_keys_locked(handle, keys):
    handle.seek(0)
    handle.truncate(0)
    for key in sorted(keys):
        handle.write(f"{key}\n")
    handle.flush()
    os.fsync(handle.fileno())


def _release_key(lock_list_path, key):
    lock_list_path.parent.mkdir(parents=True, exist_ok=True)
    lock_list_path.touch(exist_ok=True)
    with open(lock_list_path, "a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        lock_handle.seek(0)
        keys = {line.strip() for line in lock_handle if line.strip()}
        keys.discard(key)
        _write_keys_locked(lock_handle, keys)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _cleanup_old_partials(staged_dir, part_ttl_min):
    cutoff = time.time() - (part_ttl_min * 60)
    for part_path in staged_dir.glob("*.part.*"):
        try:
            if part_path.stat().st_mtime < cutoff:
                part_path.unlink(missing_ok=True)
        except FileNotFoundError:
            continue


def _rebuild_keys(staged_dir, log_path):
    keys = set()
    for r1_path in staged_dir.glob("*.R1.fastq.gz"):
        key = r1_path.name[: -len(".R1.fastq.gz")]
        r2_path = staged_dir / f"{key}.R2.fastq.gz"
        try:
            valid = (
                r1_path.exists()
                and r2_path.exists()
                and r1_path.stat().st_size > 0
                and r2_path.stat().st_size > 0
            )
        except FileNotFoundError:
            valid = False

        if valid:
            keys.add(key)
            continue

        r1_path.unlink(missing_ok=True)
        r2_path.unlink(missing_ok=True)
        _log(log_path, f"Removed incomplete staged FASTQ for {key}")

    return keys


def main():
    log_path = str(snakemake.log[0])
    src_r1 = Path(_as_single(snakemake.input.fr, "input.fr"))
    src_r2 = Path(_as_single(snakemake.input.rr, "input.rr"))
    out_r1 = Path(str(snakemake.output.fr))
    out_r2 = Path(str(snakemake.output.rr))

    lock_dir = Path(snakemake.params.get("lock_dir", ".staging_locks"))
    lock_list = lock_dir / "staged_fastq_locks.txt"
    key = f"{snakemake.wildcards.sample}.{snakemake.wildcards.lane}"
    max_staged = int(snakemake.params.get("max_staged", 6))
    wait_sec = int(snakemake.params.get("wait_sec", 30))
    part_ttl_min = int(snakemake.params.get("part_ttl_min", 120))

    staged_dir = out_r1.parent
    lock_dir.mkdir(parents=True, exist_ok=True)
    staged_dir.mkdir(parents=True, exist_ok=True)
    lock_list.touch(exist_ok=True)

    _cleanup_old_partials(staged_dir, part_ttl_min)

    with open(lock_list, "a+", encoding="utf-8") as lock_handle:
        while True:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
            keys = _rebuild_keys(staged_dir, log_path)
            if len(keys) < max_staged or key in keys:
                keys.add(key)
                _write_keys_locked(lock_handle, keys)
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
                break

            _write_keys_locked(lock_handle, keys)
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
            time.sleep(wait_sec)

    tmp_r1 = Path(f"{out_r1}.part.{os.getpid()}")
    tmp_r2 = Path(f"{out_r2}.part.{os.getpid()}")

    try:
        shutil.copyfile(src_r1, tmp_r1)
        shutil.copyfile(src_r2, tmp_r2)
        os.replace(tmp_r1, out_r1)
        os.replace(tmp_r2, out_r2)
    except Exception:
        tmp_r1.unlink(missing_ok=True)
        tmp_r2.unlink(missing_ok=True)
        _release_key(lock_list, key)
        raise

    try:
        size_bytes = out_r1.stat().st_size + out_r2.stat().st_size
        _log(log_path, f"Staged {size_bytes} bytes to local disk")
    except FileNotFoundError:
        _log(log_path, "Staging completed, but output size could not be measured")


main()
