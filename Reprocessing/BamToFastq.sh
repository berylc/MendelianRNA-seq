#!/bin/bash
picard_jar=$1
bam=$2
id=`basename $bam |sed 's/.bam//' ` 
mkdir $id
java -jar -Xmx16g $picard_jar SamToFastq \
INPUT=$bam \
OUTPUT_PER_RG=TRUE \
OUTPUT_DIR=./$id \
NON_PF=True

echo "Started gzipping"
gzip ./$id/*
echo "Finished gzipping"
