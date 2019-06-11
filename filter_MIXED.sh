#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default

source /home/wrowell/envs/Broad.env

# required: $REF $VCF
# $REF -> reference fasta; GATK sequence dict required
# $VCF -> GATK4 HaplotypeCaller vcf

echo "Extracting and filtering mixed variants from variant call file ${VCF}."
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${VCF}" \
    --output "${VCF%.*}.raw_mixed.vcf" \
    --select-type-to-include MIXED
$GATK VariantFiltration \
    --reference "${REF}" \
    --variant "${VCF%.*}.raw_mixed.vcf" \
    --output "${VCF%.*}.filtered_mixed.vcf" \
    --filter-name "DPlt3" \
    --filter-expression "DP < 3" \
    --filter-name "Qlt300" \
    --filter-expression "QUAL < 300" \
    --filter-name "QDlt2" \
    --filter-expression "QD < 2.0" \
    --filter-name "ReadPosRankSumlt-20" \
    --filter-expression "ReadPosRankSum < -20.0"
    --mask "${HPBED}" \
    --mask-extension 1 \
    --mask-name "HPmask"