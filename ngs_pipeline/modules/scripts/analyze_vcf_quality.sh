#!/bin/bash

# Script to analyze VCF quality metrics
# Usage: bash analyze_vcf_quality_fixed.sh <vcf_file>

# Activate conda environment
source /home/students/anaconda3/etc/profile.d/conda.sh
conda activate smk2

VCF=$1

if [ -z "$VCF" ]; then
    echo "Usage: $0 <vcf_file.vcf.gz>"
    exit 1
fi

BASENAME=$(basename $VCF .vcf.gz)
OUTPUT_DIR="vcf_qc_${BASENAME}"
mkdir -p $OUTPUT_DIR

echo "====================================="
echo "Analyzing: $VCF"
echo "Output directory: $OUTPUT_DIR"
echo "====================================="

# 1. Basic VCF statistics
echo -e "\n1. Basic VCF Statistics:"
bcftools stats $VCF > ${OUTPUT_DIR}/${BASENAME}_stats.txt
echo "   Saved to: ${OUTPUT_DIR}/${BASENAME}_stats.txt"

# 2. Variant counts by type
echo -e "\n2. Variant Counts by Type:"
bcftools query -f '%TYPE\n' $VCF | sort | uniq -c > ${OUTPUT_DIR}/${BASENAME}_variant_types.txt
cat ${OUTPUT_DIR}/${BASENAME}_variant_types.txt

# 3. Filter statistics
echo -e "\n3. Filter Statistics:"
bcftools query -f '%FILTER\n' $VCF | sort | uniq -c > ${OUTPUT_DIR}/${BASENAME}_filters.txt
cat ${OUTPUT_DIR}/${BASENAME}_filters.txt

# 4. Quality score distribution
echo -e "\n4. Quality Score Distribution:"
bcftools query -f '%QUAL\n' $VCF | \
    awk '{
        if ($1 < 10) q10++; 
        else if ($1 < 20) q20++; 
        else if ($1 < 30) q30++; 
        else if ($1 < 50) q50++; 
        else if ($1 < 100) q100++; 
        else q100plus++;
        total++;
    } END {
        print "QUAL < 10:    " q10 " (" int(q10/total*100) "%)";
        print "QUAL 10-20:   " q20 " (" int(q20/total*100) "%)";
        print "QUAL 20-30:   " q30 " (" int(q30/total*100) "%)";
        print "QUAL 30-50:   " q50 " (" int(q50/total*100) "%)";
        print "QUAL 50-100:  " q100 " (" int(q100/total*100) "%)";
        print "QUAL >= 100:  " q100plus " (" int(q100plus/total*100) "%)";
        print "Total:        " total;
    }' | tee ${OUTPUT_DIR}/${BASENAME}_qual_distribution.txt

# 5. Depth (DP) statistics
echo -e "\n5. Depth (DP) Statistics:"
bcftools query -f '[%DP\n]' $VCF | \
    awk '{
        if ($1 <= 5) dp5++; 
        else if ($1 <= 10) dp10++; 
        else if ($1 <= 20) dp20++; 
        else if ($1 <= 30) dp30++; 
        else if ($1 <= 50) dp50++; 
        else dp50plus++;
        sum += $1;
        count++;
    } END {
        print "DP <= 5:      " dp5 " (" int(dp5/count*100) "%)";
        print "DP 5-10:      " dp10 " (" int(dp10/count*100) "%)";
        print "DP 10-20:     " dp20 " (" int(dp20/count*100) "%)";
        print "DP 20-30:     " dp30 " (" int(dp30/count*100) "%)";
        print "DP 30-50:     " dp50 " (" int(dp50/count*100) "%)";
        print "DP > 50:      " dp50plus " (" int(dp50plus/count*100) "%)";
        print "Mean DP:      " sum/count;
        print "Total:        " count;
    }' | tee ${OUTPUT_DIR}/${BASENAME}_dp_distribution.txt

# 6. Genotype Quality (GQ) statistics
echo -e "\n6. Genotype Quality (GQ) Statistics:"
bcftools query -f '[%GQ\n]' $VCF | \
    awk '{
        if ($1 <= 10) gq10++; 
        else if ($1 <= 20) gq20++; 
        else if ($1 <= 30) gq30++; 
        else if ($1 <= 50) gq50++; 
        else gq50plus++;
        sum += $1;
        count++;
    } END {
        print "GQ <= 10:     " gq10 " (" int(gq10/count*100) "%)";
        print "GQ 10-20:     " gq20 " (" int(gq20/count*100) "%)";
        print "GQ 20-30:     " gq30 " (" int(gq30/count*100) "%)";
        print "GQ 30-50:     " gq50 " (" int(gq50/count*100) "%)";
        print "GQ > 50:      " gq50plus " (" int(gq50plus/count*100) "%)";
        print "Mean GQ:      " sum/count;
        print "Total:        " count;
    }' | tee ${OUTPUT_DIR}/${BASENAME}_gq_distribution.txt

# 7. Genotype distribution
echo -e "\n7. Genotype Distribution:"
bcftools query -f '[%GT\n]' $VCF | sort | uniq -c | \
    awk '{print $2 ": " $1}' | tee ${OUTPUT_DIR}/${BASENAME}_genotypes.txt

# 8. Hard filter metrics (QD, FS, MQ, etc.)
echo -e "\n8. Hard Filter Metrics:"

# QD (Quality by Depth)
echo "QD (Quality by Depth):"
bcftools query -f '%INFO/QD\n' $VCF 2>/dev/null | \
    awk '{
        if ($1 < 2) qd2++;
        sum += $1;
        count++;
    } END {
        if (count > 0) {
            print "  QD < 2:     " qd2 " (" int(qd2/count*100) "%)";
            print "  Mean QD:    " sum/count;
        } else {
            print "  QD not available";
        }
    }' | tee -a ${OUTPUT_DIR}/${BASENAME}_hard_filters.txt

# FS (Fisher Strand)
echo "FS (Fisher Strand):"
bcftools query -f '%INFO/FS\n' $VCF 2>/dev/null | \
    awk '{
        if ($1 > 60) fs60++;
        sum += $1;
        count++;
    } END {
        if (count > 0) {
            print "  FS > 60:    " fs60 " (" int(fs60/count*100) "%)";
            print "  Mean FS:    " sum/count;
        } else {
            print "  FS not available";
        }
    }' | tee -a ${OUTPUT_DIR}/${BASENAME}_hard_filters.txt

# MQ (Mapping Quality)
echo "MQ (Mapping Quality):"
bcftools query -f '%INFO/MQ\n' $VCF 2>/dev/null | \
    awk '{
        if ($1 < 40) mq40++;
        sum += $1;
        count++;
    } END {
        if (count > 0) {
            print "  MQ < 40:    " mq40 " (" int(mq40/count*100) "%)";
            print "  Mean MQ:    " sum/count;
        } else {
            print "  MQ not available";
        }
    }' | tee -a ${OUTPUT_DIR}/${BASENAME}_hard_filters.txt

# 9. Transition/Transversion ratio
echo -e "\n9. Ti/Tv Ratio:"
bcftools stats $VCF | grep "TSTV" | head -1 | \
    awk '{print "  Ti/Tv ratio: " $5}' | tee ${OUTPUT_DIR}/${BASENAME}_titv.txt

# 10. Chromosomal distribution
echo -e "\n10. Chromosomal Distribution:"
bcftools query -f '%CHROM\n' $VCF | sort | uniq -c | \
    sort -rn | head -25 > ${OUTPUT_DIR}/${BASENAME}_chr_distribution.txt
cat ${OUTPUT_DIR}/${BASENAME}_chr_distribution.txt

echo -e "\n====================================="
echo "Analysis complete!"
echo "All results saved to: $OUTPUT_DIR/"
echo "====================================="

# Create summary report
cat > ${OUTPUT_DIR}/${BASENAME}_summary.txt << EOF
VCF Quality Analysis Summary
============================
File: $VCF
Date: $(date)

Key Metrics:
-----------
$(grep "number of records:" ${OUTPUT_DIR}/${BASENAME}_stats.txt | head -1)
$(grep "number of SNPs:" ${OUTPUT_DIR}/${BASENAME}_stats.txt | head -1)
$(grep "number of indels:" ${OUTPUT_DIR}/${BASENAME}_stats.txt | head -1)

Filter Status:
-------------
$(cat ${OUTPUT_DIR}/${BASENAME}_filters.txt)

Quality Distribution:
--------------------
$(cat ${OUTPUT_DIR}/${BASENAME}_qual_distribution.txt)

Depth Distribution:
------------------
$(cat ${OUTPUT_DIR}/${BASENAME}_dp_distribution.txt)

GQ Distribution:
---------------
$(cat ${OUTPUT_DIR}/${BASENAME}_gq_distribution.txt)

Ti/Tv Ratio:
-----------
$(cat ${OUTPUT_DIR}/${BASENAME}_titv.txt)
EOF

echo -e "\nSummary report: ${OUTPUT_DIR}/${BASENAME}_summary.txt"   