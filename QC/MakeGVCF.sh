#!/bin/bash
gatk_jar=$1
hg_fasta=$2
bam=$3

id=`basename $bam | sed 's/.bam//'`
java -Xmx4000m -jar $gatk_jar \
-T HaplotypeCaller --disable_auto_index_creation_and_locking_when_reading_rods \
-ERC GVCF --max_alternate_alleles 3 -variant_index_parameter 128000 -variant_index_type LINEAR --read_filter OverclippedRead \
-contamination 0.0 \
-R $hg_fasta \
-U ALLOW_N_CIGAR_READS \
-L ../data/purcell.interval_list \
-I $bam \
-o ${id}.vcf.gz
