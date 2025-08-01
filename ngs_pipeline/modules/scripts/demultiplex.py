import os
import subprocess
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import logging
import shutil

# Set up logging
logging.basicConfig(level=logging.INFO,
                   format='%(asctime)s - %(levelname)s - %(message)s')

def find_bcl_directories(root_dir):
    """
    Recursively find all directories containing BCL files
    
    Parameters:
        root_dir (str): Root directory to start search from
        
    Returns:
        list: List of directories containing BCL files
    """
    root_path = Path(root_dir)
    bcl_directories = set(bcl_file.parent for bcl_file in root_path.rglob("*.bcl"))
    logging.info(f"Found {len(bcl_directories)} BCL directories")
    return list(bcl_directories)

def convert_single_bcl(bcl_dir, sample_sheet, output_base_dir, threads=None):
    """
    Convert BCL files from a directory to FASTQ
    
    Parameters:
        bcl_dir (Path): Directory containing BCL files
        sample_sheet (Path): Path to SampleSheet.csv
        output_base_dir (str): Base directory for output
        threads (int, optional): Number of threads for processing
    """
    # Create output directory matching input structure
    relative_path = bcl_dir.name
    output_dir = Path(output_base_dir) / relative_path
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Build command
    cmd = ['bcl2fastq',
           '--input-dir', str(bcl_dir),
           '--output-dir', str(output_dir),
           '--sample-sheet', str(sample_sheet),
           '--no-lane-splitting']
    
    if threads:
        cmd.extend(['--processing-threads', str(threads)])
    
    try:
        logging.info(f"Processing {bcl_dir}")
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logging.info(f"Successfully converted {bcl_dir}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error converting {bcl_dir}: {e.stderr}")
        return False

def batch_convert_bcl(input_dir, output_dir, sample_sheet, max_workers=4, threads_per_job=1):
    """
    Convert all BCL files found in the directory structure using a single SampleSheet
    
    Parameters:
        input_dir (str): Root input directory
        output_dir (str): Root output directory
        sample_sheet (str): Path to SampleSheet.csv
        max_workers (int): Maximum number of parallel conversions
        threads_per_job (int): Threads per conversion job
    """
    # Verify sample sheet exists
    sample_sheet_path = Path(sample_sheet)
    if not sample_sheet_path.exists():
        logging.error(f"SampleSheet not found: {sample_sheet}")
        return
    
    # Find all BCL directories
    bcl_dirs = find_bcl_directories(input_dir)
    
    if not bcl_dirs:
        logging.error("No BCL directories found!")
        return
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Convert files in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                convert_single_bcl,
                bcl_dir,
                sample_sheet_path,
                output_dir,
                threads_per_job
            )
            for bcl_dir in bcl_dirs
        ]
        
        # Wait for all conversions to complete
        successful = 0
        for future in futures:
            if future.result():
                successful += 1
    
    logging.info(f"Conversion complete. Successfully converted {successful}/{len(bcl_dirs)} directories")

def main():
    parser = argparse.ArgumentParser(description='Batch convert BCL files to FASTQ format')
    
    parser.add_argument('-i', '--input', required=True,
                      help='Root directory containing BCL files in subdirectories')
    parser.add_argument('-o', '--output', required=True,
                      help='Output directory for FASTQ files')
    parser.add_argument('-s', '--sample-sheet', required=True,
                      help='Path to SampleSheet.csv file')
    parser.add_argument('-w', '--workers', type=int, default=4,
                      help='Number of parallel conversion jobs (default: 4)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                      help='Number of threads per conversion job (default: 1)')
    
    args = parser.parse_args()
    
    batch_convert_bcl(
        input_dir=args.input,
        output_dir=args.output,
        sample_sheet=args.sample_sheet,
        max_workers=args.workers,
        threads_per_job=args.threads
    )

if __name__ == "__main__":
    main()
