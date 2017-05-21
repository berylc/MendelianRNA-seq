#!/bin/bash
tophatBam=$1
fqDir=$2
starPath=$3
starGenomeFile=$4
junctionFile=$5
#Giving it the original TopHat bam, to extract sample name and get readgroups. 
sampname=`basename $tophatBam | sed 's/.bam//'`
fqPath=$sampname

allpair1=`ls -1 $fqDir/$fqPath/*_1.fastq.gz | sort | tr "\n" "," |  sed 's/.$//'`
allpair2=`ls -1 $fqDir/$fqPath/*_2.fastq.gz | sort | tr "\n" "," |  sed 's/.$//'`

#Only using ID from Read group
samtools view -H $tophatBam  | egrep "^@RG" | cut -f2 | sort | sed ':a;N;$!ba;s/\n/ , /g' | tr "\t" " " > ReadGroups_${sampname}.txt

#Star command
./GeneralAlignment.sh \
$allpair1 \
$allpair2 \
$starGenomeFile \
str_PE \
8 \
${sampname}_2ndPass \
76 \
- \
ReadGroups_${sampname}.txt \
yes \
- \
readgroup \
sorted \
$starPath
