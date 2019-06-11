#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 4

source /mnt/software/Modules/current/init/bash
module load samtools
source Broad.env

# required: $SAMPLE $REF $BAMLIST $BED $PCRINDELMODEL $PLOIDY $MAPQ
# $SAMPLE -> sample name
# $REF -> reference fasta; GATK sequence dict required
# $BAMLIST -> either single BAM or list of BAM file names (e.g. bam.list)
# $BED -> BED file describing targeted region
# $PCRINDELMODEL -> [NONE, CONSERVATIVE, AGGRESSIVE, HOSTILE]
#                   See https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.6.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php#--pcr-indel-model
#                   Tweak depending on chemistry and whether amplification was used.  HOSTILE works well with amplified DNA and 2.0 Chemistry.
# $PLOIDY -> [1, 2]
# $MAPQ -> for pbmm2/minimap2, [1..60]

echo "Calling variants for ${BED} reads from ${BAMLIST} against ${REF}, with pcr-indel-model ${PCRINDELMODEL}, minimum MAPQ ${MAPQ}, and ploidy expectation ${PLOIDY}."
$GATK HaplotypeCaller \
 --reference "${REF}" \
 --input "${BAMLIST}" \
 --output "${SAMPLE}.${PLOIDY}ploid.${PCRINDELMODEL}.HaplotypeCaller.g.vcf.gz" \
 --sample-ploidy "${PLOIDY}" \
 --pcr-indel-model "${PCRINDELMODEL}" \
 --intervals "${BED}" \
 --read-filter MappingQualityReadFilter \
 --read-filter NotSecondaryAlignmentReadFilter \
 --read-filter NotSupplementaryAlignmentReadFilter \
 --minimum-mapping-quality "${MAPQ}" \
 --emit-ref-confidence GVCF \
 --annotation-group StandardAnnotation \
 --annotation-group AS_StandardAnnotation \
 --annotation-group StandardHCAnnotation
