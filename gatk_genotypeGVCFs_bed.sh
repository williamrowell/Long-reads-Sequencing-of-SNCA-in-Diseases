#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 4

source /home/wrowell/envs/Broad.env

# required: $GVCF $REF $BED
# $GVCF -> *.g.vcf.gz file produced by GATK4 HaplotypeCaller using the parameters in gatk_hc_gvcf_bed.sh
# $REF -> reference fasta; GATK sequence dict required
# $BED -> BED file describing targeted region

echo "Calling variants for ${BED} reads from ${GVCF} against ${REF}."
$GATK GenotypeGVCFs \
 --reference "${REF}" \
 --dbsnp "${DBSNP}" \
 --variant "${GVCF}" \
 --output "${GVCF%.g.vcf.gz}.vcf.gz" \
 --intervals "${BED}" \
 --annotation-group StandardAnnotation \
 --annotation-group AS_StandardAnnotation \
 --annotation-group StandardHCAnnotation
