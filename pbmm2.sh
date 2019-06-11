#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 16

source /mnt/software/Modules/current/init/bash
module load pbmm2

# Produce minimap2 aligned BAMs that conform to PacBio BAM specifications.
# required: $SAMPLE $REF $CCSSET $OUTSET
# $SAMPLE -> sample name
# $REF -> reference fasta or referenceset.xml
# $CCSSET -> either *.consensusreadset.xml or *.consensusreads.bam
# $OUTSET -> either *.consensusalignmentset.xml or *.consensusreads.aln.bam

echo "Aligning ${CCSSET} to ${REF} to produce ${OUTSET}."
pbmm2 align --alignment-threads 14 \
	    --sort-threads 2 \
	    --sort-memory 16G \
	    --preset CCS \
	    --sort \
	    --sample ${SAMPLE} \
	    ${CCSSET} \
	    ${REF} \
	    ${OUTSET}
