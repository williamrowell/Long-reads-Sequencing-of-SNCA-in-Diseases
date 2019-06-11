#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default
#$ -pe smp 16

source /mnt/software/Modules/current/init/bash
module load smrtanalysis/6.0.0

# required: $REF $SVSIGFOFN
# $REF -> reference fasta
# $SVSIGFOFN -> file with one *.svsig.gz per line

refname="${REF##*/}"
svsigname="${SVSIGFOFN##*/}"

echo "Calling structural variants for svsig files in ${SVSIGFOFN} using the reference ${REF}."

pbsv call --num-threads 16 \
        "${REF}" \
        "${SVSIGFOFN}" \
        "${svsigname%.svsig.gz}.var.vcf"