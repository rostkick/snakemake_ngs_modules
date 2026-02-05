#!/usr/bin/env python3
"""
VCF File Diagnostic Tool
"""

import gzip
import sys
import pandas as pd
from collections import Counter

def parse_vcf(vcf_file):
    """Parse VCF file and collect statistics"""
    
    stats = {
        'total_variants': 0,
        'chromosomes': Counter(),
        'filters': Counter(),
        'qualities': [],
        'formats': Counter(),
        'samples': set(),
        'dp_values': [],
        'gq_values': [],
        'genotypes': Counter(),
        'info_fields': Counter()
    }
    
    print(f"Analyzing VCF file: {vcf_file}")
    print("="*60)
    
    # Determine if file is gzipped
    if vcf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    
    with opener(vcf_file, 'rt') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and headers
            if not line or line.startswith('##'):
                continue
            
            # Parse header line
            if line.startswith('#CHROM'):
                header = line[1:].split('\t')  # Remove # from #CHROM
                stats['samples'] = header[9:]  # Samples start at column 10
                print(f"Samples found: {len(stats['samples'])}")
                if stats['samples']:
                    print(f"Sample names: {', '.join(stats['samples'])}")
                continue
            
            # Parse variant line
            stats['total_variants'] += 1
            fields = line.split('\t')
            
            if len(fields) < 8:
                continue
            
            # Basic fields
            chrom = fields[0]
            pos = fields[1]
            var_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            qual = fields[5]
            vcf_filter = fields[6]
            info = fields[7]
            format_col = fields[8] if len(fields) > 8 else ""
            sample_data = fields[9:] if len(fields) > 9 else []
            
            # Collect statistics
            stats['chromosomes'][chrom] += 1
            
            # Filter field
            if vcf_filter == '.':
                vcf_filter = 'PASS'
            stats['filters'][vcf_filter] += 1
            
            # Quality
            if qual != '.':
                try:
                    stats['qualities'].append(float(qual))
                except:
                    pass
            
            # Format field
            stats['formats'][format_col] += 1
            
            # INFO field analysis
            if info != '.':
                info_fields = info.split(';')
                for field in info_fields:
                    if '=' in field:
                        key = field.split('=')[0]
                        stats['info_fields'][key] += 1
                    else:
                        stats['info_fields'][field] += 1
            
            # Sample data analysis
            for sample in sample_data:
                if format_col and sample != '.':
                    format_keys = format_col.split(':')
                    sample_values = sample.split(':')
                    
                    for key, value in zip(format_keys, sample_values):
                        if key == 'GT':
                            stats['genotypes'][value] += 1
                        elif key == 'DP':
                            if value != '.':
                                try:
                                    stats['dp_values'].append(int(value))
                                except:
                                    pass
                        elif key == 'GQ':
                            if value != '.':
                                try:
                                    stats['gq_values'].append(int(value))
                                except:
                                    pass
    
    return stats

def print_statistics(stats):
    """Print collected statistics"""
    
    print(f"\nTotal variants: {stats['total_variants']:,}")
    
    # Chromosome distribution
    print("\n" + "="*60)
    print("CHROMOSOME DISTRIBUTION")
    print("="*60)
    
    total_chrom = sum(stats['chromosomes'].values())
    for chrom, count in sorted(stats['chromosomes'].items()):
        percentage = count / total_chrom * 100
        print(f"{chrom:10}: {count:8,} ({percentage:5.1f}%)")
    
    # Filter distribution
    print("\n" + "="*60)
    print("FILTER FIELD DISTRIBUTION")
    print("="*60)
    
    for filt, count in stats['filters'].most_common():
        percentage = count / stats['total_variants'] * 100
        print(f"{filt:20}: {count:8,} ({percentage:5.1f}%)")
    
    # Quality statistics
    if stats['qualities']:
        print("\n" + "="*60)
        print("QUAL (QUALITY) FIELD STATISTICS")
        print("="*60)
        
        qual_series = pd.Series(stats['qualities'])
        print(f"Mean QUAL: {qual_series.mean():.1f}")
        print(f"Median QUAL: {qual_series.median():.1f}")
        print(f"Min QUAL: {qual_series.min():.1f}")
        print(f"Max QUAL: {qual_series.max():.1f}")
        
        # Quality thresholds
        thresholds = [1, 10, 20, 30, 50, 100, 200]
        print(f"\nQUAL Thresholds:")
        for threshold in thresholds:
            below = (qual_series < threshold).sum()
            percentage = below / len(qual_series) * 100
            print(f"  QUAL < {threshold:3}: {below:8,} ({percentage:5.1f}%)")
    
    # Format field analysis
    print("\n" + "="*60)
    print("FORMAT FIELD ANALYSIS")
    print("="*60)
    
    for fmt, count in stats['formats'].most_common(10):  # Top 10 formats
        percentage = count / stats['total_variants'] * 100
        print(f"{fmt:40}: {count:8,} ({percentage:5.1f}%)")
    
    # Genotype distribution
    if stats['genotypes']:
        print("\n" + "="*60)
        print("GENOTYPE (GT) DISTRIBUTION")
        print("="*60)
        
        total_geno = sum(stats['genotypes'].values())
        for gt, count in stats['genotypes'].most_common():
            percentage = count / total_geno * 100
            print(f"{gt:10}: {count:8,} ({percentage:5.1f}%)")
    
    # DP statistics
    if stats['dp_values']:
        print("\n" + "="*60)
        print("DEPTH (DP) STATISTICS")
        print("="*60)
        
        dp_series = pd.Series(stats['dp_values'])
        print(f"Mean DP: {dp_series.mean():.1f}")
        print(f"Median DP: {dp_series.median():.1f}")
        print(f"Min DP: {dp_series.min():.1f}")
        print(f"Max DP: {dp_series.max():.1f}")
        
        # DP distribution
        print(f"\nDP Distribution:")
        thresholds = [1, 5, 10, 20, 30, 50, 100]
        for threshold in thresholds:
            below = (dp_series < threshold).sum()
            at_or_below = (dp_series <= threshold).sum()
            pct_below = below / len(dp_series) * 100
            pct_at_or_below = at_or_below / len(dp_series) * 100
            print(f"  DP < {threshold:3}: {below:8,} ({pct_below:5.1f}%)")
            print(f"  DP ≤ {threshold:3}: {at_or_below:8,} ({pct_at_or_below:5.1f}%)")
    
    # GQ statistics
    if stats['gq_values']:
        print("\n" + "="*60)
        print("GENOTYPE QUALITY (GQ) STATISTICS")
        print("="*60)
        
        gq_series = pd.Series(stats['gq_values'])
        print(f"Mean GQ: {gq_series.mean():.1f}")
        print(f"Median GQ: {gq_series.median():.1f}")
        print(f"Min GQ: {gq_series.min():.1f}")
        print(f"Max GQ: {gq_series.max():.1f}")
        
        # GQ distribution
        print(f"\nGQ Distribution:")
        thresholds = [1, 5, 10, 20, 30, 50, 99]
        for threshold in thresholds:
            below = (gq_series < threshold).sum()
            at_or_below = (gq_series <= threshold).sum()
            pct_below = below / len(gq_series) * 100
            pct_at_or_below = at_or_below / len(gq_series) * 100
            print(f"  GQ < {threshold:3}: {below:8,} ({pct_below:5.1f}%)")
            print(f"  GQ ≤ {threshold:3}: {at_or_below:8,} ({pct_at_or_below:5.1f}%)")
    
    # INFO fields
    if stats['info_fields']:
        print("\n" + "="*60)
        print("TOP 20 INFO FIELDS")
        print("="*60)
        
        for field, count in stats['info_fields'].most_common(20):
            percentage = count / stats['total_variants'] * 100
            print(f"{field:30}: {count:8,} ({percentage:5.1f}%)")

def main():
    if len(sys.argv) != 2:
        print("Usage: python vcf_diagnostic.py input.vcf[.gz]")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    
    try:
        stats = parse_vcf(vcf_file)
        print_statistics(stats)
    except Exception as e:
        print(f"Error analyzing VCF file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()