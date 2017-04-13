#!/bin/bash
loc=$1
id=`basename $loc |sed 's/.bam//' ` 
mkdir $id
java -jar -Xmx16g /path/to/picard.jar SamToFastq \
INPUT=${loc} \
OUTPUT_PER_RG=TRUE \
OUTPUT_DIR=./$id \
NON_PF=True

echo "Started gzipping"
gzip ./$id/*
echo "Finished gzipping"
