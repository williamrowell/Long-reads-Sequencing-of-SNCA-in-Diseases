#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 16

source /mnt/software/Modules/current/init/bash
module load smrtanalysis/6.0.0 minimap2 samtools

# Align one read per ZMW to reference using minimap2
# required: $MOVIE $REF $SAMPLE
# $MOVIE -> subreadset.xml
# $REF -> reference fasta
# $SAMPLE -> sample name

moviename="${MOVIE##*/}"
refname="${REF##*/}"
RGRECORD="@RG\tID:${SAMPLE}_${moviename%.*}\tSM:${SAMPLE}"

echo "Aligning ${MOVIE} to ${REF}."
echo "Adding @RG record with ID:${SAMPLE}_${moviename%.*} and SM:${SAMPLE}"

pbsv fasta ${MOVIE} | \
        minimap2 -t 16 \
                -ax map-pb \
                --eqx -L \
                -O 5,56 \
                -E 4,1 \
                -B 5 \
                --secondary=no \
                -z 400,50 \
                -r 2k -Y \
                -R "${RGRECORD}" \
                "${REF}" - | \
        samtools sort --threads 2 - > "${refname%.*}.${moviename%.*}.aln.bam"
samtools index -@ 15 "${refname%.*}.${moviename%.*}.aln.bam"