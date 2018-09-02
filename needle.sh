#!/bin/bash
#========
# USAGE: count_one_fraction.sh reference.fa output_name activity
#========

referenceName="$1"
referenceFileBase=${1%.fa}
referenceSequence=$(sed -n '2p' $1)
baseName="$2"
activity="$3"

set -e

## next need to do alignment with needle-all
#if [ -s $baseName.$activity.aln ]
#then
#    rm $baseName.$activity.aln
#fi
#
#needleall -gapopen 15 -gapextend 0.5 -asequence $referenceFileBase.fa -supper1 \
#-bsequence $baseName.$activity.fa -supper2 \
#-aformat3 fasta -outfile $baseName.$activity.aln -errfile $baseName.$activity.err
#echo -e "Multiple to one alignment complete\n"

echo -e "Done with set $baseName.$activity"
