import argparse
import fcntl
from pathlib import Path


def write_keys_locked(handle, keys):
    handle.seek(0)
    handle.truncate(0)
    for key in sorted(keys):
        handle.write(f"{key}\n")
    handle.flush()


def release_key(lock_dir, key, log_path=None):
    lock_dir_path = Path(lock_dir)
    lock_list = lock_dir_path / "staged_fastq_locks.txt"
    lock_dir_path.mkdir(parents=True, exist_ok=True)
    lock_list.touch(exist_ok=True)

    with open(lock_list, "a+", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        keys = {line.strip() for line in handle if line.strip()}
        keys.discard(key)
        write_keys_locked(handle, keys)
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)

    if log_path:
        with open(log_path, "a", encoding="utf-8") as log_handle:
            log_handle.write(f"Released staging lock for {key}\n")


def parse_args():
    parser = argparse.ArgumentParser(description="Release staged FASTQ lock entry")
    parser.add_argument("--lock-dir", default=".staging_locks")
    parser.add_argument("--key", required=True)
    parser.add_argument("--log", default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    release_key(args.lock_dir, args.key, args.log)


if __name__ == "__main__":
    main()
