#!/bin/bash
bam=$1
id=`basename $bam | sed 's/.bam//'`
java -Xmx4000m -jar "/path/to/jar/gatk.jar" \
-T HaplotypeCaller --disable_auto_index_creation_and_locking_when_reading_rods \
-ERC GVCF --max_alternate_alleles 3 -variant_index_parameter 128000 -variant_index_type LINEAR --read_filter OverclippedRead \
-contamination 0.0 \
-R "/path/to/Homo_sapiens_assembly19.fasta" \
-U ALLOW_N_CIGAR_READS \
-L ../purcell.interval_list \
-I $bam \
-o ${id}.vcf.gz
