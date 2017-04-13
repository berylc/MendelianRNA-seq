#!/bin/bash
gvcflist=$1
java -Xmx28g -jar /seq/software/picard/current/bin/GenomeAnalysisTK-3.4-g3c929b0.jar \
-T GenotypeGVCFs --disable_auto_index_creation_and_locking_when_reading_rods \
-R /seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta \
-L ../data/purcell.interval_list \
-V $gvcflist \
-o out.joint.vcf.gz
