#!/bin/bash
loc=$1
fqDir=$2
#Giving it the original bam name, just to extract sample name.
sampname=`basename $loc | sed 's/.bam//'`
fqPath=./$sampname

#FqPath is just sampname in the folder
allpair1=`ls -1 $fqDir$fqPath/*_1.fastq.gz | sort | tr "\n" "," |  sed 's/.$//'`
allpair2=`ls -1 $fqDir$fqPath/*_2.fastq.gz | sort | tr "\n" "," |  sed 's/.$//'`

#Only using ID from Read group
samtools view -H $loc  | egrep "^@RG" | cut -f2 | sort | sed ':a;N;$!ba;s/\n/ , /g' | tr "\t" " " > ReadGroups_${sampname}.txt

#Star command
./GeneralAlignment.sh \
$allpair1 \
$allpair2 \
/path/to/genome/file/ \
str_PE \
8 \
${sampname}_1stPass \
101 \
- \
ReadGroups_${sampname}.txt \
yes \
- \
readgroup \
sorted

