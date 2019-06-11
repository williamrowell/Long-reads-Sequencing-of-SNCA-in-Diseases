#!/bin/bash
#$ -cwd
#$ -V
#$ -j y
#$ -S /bin/bash
#$ -q default

source /mnt/software/Modules/current/init/bash
module load smrtanalysis/6.0.0

echo "Discovering structural variant signatures in ${BAM}."
# required: $BAM
# $BAM -> BAM produced by aligning a single subread from each ZMW to a reference with minimap2

pbsv discover "${BAM}" "${BAM%.*}.svsig.gz"