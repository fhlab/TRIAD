#!/bin/bash
#========
# USAGE: count.sh reference.fa output_name activity
#========

referenceName="$1"
referenceFileBase=${1%.fa}
referenceSequence=$(sed -n '2p' $1)
baseName="$2"
activity="$3"


set -e

# create two SAM files and extract properly mapped reads that don't contain N
echo -e "\nMapping valid reads to reference..."
# bowtie2-build $referenceName $referenceFileBase
echo -e "Mapping unassembled reads\n"
/opt/bowtie2-2.3.0/bowtie2 -x $referenceFileBase -1 $baseName.$activity.unassembled.forward.fastq \
-2 $baseName.$activity.unassembled.reverse.fastq -S $baseName.$activity.unassembled.sam
echo -e "\nMapping assembled reads"
/opt/bowtie2-2.3.0/bowtie2 -x $referenceFileBase -U $baseName.$activity.assembled.fastq -S $baseName.$activity.assembled.sam


# sort the sam files to get depth per position
echo -e "Calculating coverage per position..."
samtools sort -o $baseName.$activity.unassembled.sorted.sam $baseName.$activity.unassembled.sam
samtools depth -a -m 1000000 $baseName.$activity.unassembled.sorted.sam > $baseName.$activity.unassembled.depth.txt
samtools sort -o $baseName.$activity.assembled.sorted.sam $baseName.$activity.assembled.sam
samtools depth -a -m 1000000 $baseName.$activity.assembled.sorted.sam > $baseName.$activity.assembled.depth.txt
echo -e "Bowtie2 alignment complete"

 Takes properly aligned matches and extracts reads
 I want columns 1 (name) and 10 (read)

echo -e "\nExtracting reads from SAM files..."
if [ -s $baseName.$activity.fa ]
then
    rm $baseName.$activity.fa
fi
cat $baseName.$activity.assembled.sam | grep -E "^\S+[	](0|16)" | cut -f 1,10 | \
awk 'BEGIN {OFS="\n"} {if ($2~/[AGCTN]{20,}/) {print "\n>"$1,$2}}' > $baseName.$activity.blank.fa
cat $baseName.$activity.unassembled.sam | grep -E "^\S+[	](99|147|83|163)" | cut -f 1,10 | \
awk 'BEGIN {OFS="\n"} {if ($2~/[AGCTN]{20,}/) {print "\n>"$1,$2}}' >> $baseName.$activity.blank.fa
awk 'NF' $baseName.$activity.blank.fa > $baseName.$activity.fa
echo -e "Extracting proper reads complete"

# next need to do alignment with needle-all
if [ -s $baseName.$activity.aln ]
then
    rm $baseName.$activity.aln
fi

needleall -gapopen 15 -gapextend 0.5 -asequence $referenceFileBase.fa -supper1 \
-bsequence $baseName.$activity.fa -supper2 \
-aformat3 fasta -outfile $baseName.$activity.aln -errfile $baseName.$activity.err
echo -e "Multiple to one alignment complete\n"

echo -e "Done with set $baseName.$activity"
