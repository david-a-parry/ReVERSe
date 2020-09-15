#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
VCF=$DIR/test_data/test_in.vcf.gz
PED=$DIR/test_data/test.ped
OUTDIR=$DIR/test_output
mkdir -p $OUTDIR
PREFIX=${OUTDIR}/test_output

set -euo pipefail

echo $(date) Running ReVERSe_count
ReVERSe_count -i $VCF \
    --ped $PED \
    -e $DIR/test_data/test_gnomad_wes.vcf.gz \
    -g $DIR/test_data/test_gnomad_wgs.vcf.gz \
    --exome_coverage_files $DIR/test_data/test_wes_cov.tsv.gz \
    --genome_coverage_files $DIR/test_data/test_wgs_cov.tsv.gz \
    --gnomad_version 2.1 \
    -v $PREFIX.rev_counts.vcf.gz \
    -t $PREFIX.rev_counts.txt.gz 

echo $(date) Running ReVERSe_seg
ReVERSe_seg -i $PREFIX.rev_counts.vcf.gz \
    --ped $PED \
    --freq 0.005 \
    --canonical \
    --pops afr amr eas fin nfe sas \
    -v 1e-5 \
    --impact HIGH \
    -o $PREFIX.rev_seg.vcf.gz

echo $(date) Running ReVERSe_reporter
ReVERSe_reporter $PREFIX.rev_seg.vcf.gz \
    $PREFIX.rev_seg.report.csv \
    $PED

echo $(date) Finished running tests - output in $OUTDIR
echo $(date) Exit status:
echo $?

