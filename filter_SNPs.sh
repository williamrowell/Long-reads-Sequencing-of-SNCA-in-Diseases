#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default

source Broad.env

# required: $REF $VCF
# $REF -> reference fasta; GATK sequence dict required
# $VCF -> GATK4 HaplotypeCaller vcf

echo "Extracting and filtering SNPs from variant call file ${VCF}."
$GATK SelectVariants \
    --reference "${REF}" \
    --variant "${VCF}" \
    --output "${VCF%.*}.raw_snps.vcf" \
    --select-type-to-include SNP
$GATK VariantFiltration \
    --reference "${REF}" \
    --variant "${VCF%.*}.raw_snps.vcf" \
    --output "${VCF%.*}.filtered_snps.vcf" \
    --filter-expression "DP < 3" \
    --filter-name "DPlt3" \
    --filter-expression "QUAL < 300" \
    --filter-name "Qlt300" \
    --filter-expression "QD < 2.0" \
    --filter-name "QDlt2" \
    --filter-expression "MQ < 30.0" \
    --filter-name "MQlt30" \
    --filter-expression "MQRankSum < -12.5" \
    --filter-name "MQRankSumlt-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" \
    --filter-name "ReadPosRankSumlt-8"