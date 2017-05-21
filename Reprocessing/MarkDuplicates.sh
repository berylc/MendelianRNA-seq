#!/bin/bash
picard_jar=$1
bam=$2
id=`basename $bam | sed 's/_2ndPassAligned.sortedByCoord.out.bam//'`

java -jar -Xmx20g $picard_jar MarkDuplicates \
INPUT=$bam \
OUTPUT=./${id}.sorted.deduped.bam \
METRICS_FILE=./${id}.sorted.deduped.metrics \
CREATE_INDEX=TRUE \
SORTING_COLLECTION_SIZE_RATIO=0.1 \
ASSUME_SORTED=TRUE \
VALIDATION_STRINGENCY=LENIENT
