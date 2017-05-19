
#!/bin/bash

# STAR mapping pipeline (can be used for both 1-pass and all-sample 2-pass runs by including the junctions from 1-pass run)
# usage: from an empty working directory, run
# ./STAR.sh (read1) (read2 or "") (STARgenomeDir) (dataType) (nThreadsSTAR) (sampleName) (readLength) (annotationFile) (readGroupAndOtherCommentsFile) (readFilesCommand) (junctionsFrom1Pass)

read1=$1 
#fastq file for read1
read2=$2 
#fastq file for read1, use "" if single-end
STARgenomeDir=$3 
dataType=$4 
#RNA-seq type, possible values: str_SE str_PE unstr_SE unstr_PE
nThreadsSTAR=$5 
#number of threads for STAR
sname=$6
#output prefix
readLength=$7 
#possible values 25,50,76,125,150
annotation=$8 
#annotation file, use "-" if no annotation, annotations not needed if used at the genome generation step
readGroups=$9
#file of read groups
zipped=${10} 
# yes for gzipped fastq files, no for unzipped
j=${11} 
# path(s) to junction file(s)
fqType=${12} 
#How are your fqs split up? possible values: readgroup or pairs
outBamType=${13} 
#output BAM format. possible values: sorted or unsorted
starExecutable=${14}
#path to STAR, e.g. STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR


case $j in
-)
	junctions=${11}
	save="Basic"
	;;
*)
	junctions=${11} # use "-" if 1-pass run
	save="All"
	;;
esac

STAR=$starExecutable

# STAR parameters: common for all runs
STARparCommon=" --genomeDir $STARgenomeDir --readFilesIn $read1 $read2 --outSAMunmapped Within --outSAMmapqUnique 60 \
	--outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFileNamePrefix ./$sname \
	--alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --sjdbGTFfile $annotation \
	--sjdbFileChrStartEnd $junctions --limitSjdbInsertNsj 3000000 --sjdbInsertSave $save"

# STAR parameters: by read length
case "$readLength" in
76)
	STARparReadLength=" --outFilterMismatchNoverReadLmax 0.1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 \
	--sjdbOverhang 75 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --alignSoftClipAtReferenceEnds No --chimJunctionOverhangMin 15"
	;;
50) 	STARparReadLength=" --outFilterMismatchNoverReadLmax 0.08 --alignSJoverhangMin 6 --alignSJDBoverhangMin 1 --sjdbScore 1  \
	--sjdbOverhang 49"
	;;
25)
	STARparReadLength=" --outFilterMismatchNoverReadLmax 0.16 --alignSJoverhangMin 6 --alignSJDBoverhangMin 3 --sjdbScore 3 \
	--sjdbOverhang 24 --outSJFilterOverhangMin 10 8 8 8"
	;;
101)
	STARparReadLength=" --outFilterMismatchNoverReadLmax 0.07 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 \
	--sjdbOverhang 100 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --alignSoftClipAtReferenceEnds No --chimJunctionOverhangMin 15"
	;;
	
125)
	STARparReadLength=" --outFilterMismatchNoverReadLmax 0.1 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 \
	--sjdbOverhang 124 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --alignSoftClipAtReferenceEnds No --chimJunctionOverhangMin 15"
	;;

150)
	STARparReadLength=" --outFilterMismatchNoverReadLmax 0.05 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1  \
	--sjdbOverhang 149 --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --alignSoftClipAtReferenceEnds No --chimJunctionOverhangMin 15"
	;;

esac

# STAR parameters: zipping of fastq files

case "$zipped" in
no)
	STARparZip=" --readFilesCommand cat"
	;;
yes)
	STARparZip=" --readFilesCommand zcat"
	;;
esac

# STAR parameters: run-time, controlled by DCC
STARparRun=" --runThreadN $nThreadsSTAR --genomeLoad NoSharedMemory --limitBAMsortRAM 100000000000"


# STAR parameters: type of BAM output: quantification or sorted BAM or both
#   OPTION: sorted BAM output
#STARparBAM="--outSAMtype BAM SortedByCoordinate"
#   OPTION: transcritomic BAM for quantification
## STARparBAM="--outSAMtype None --quantMode TranscriptomeSAM"
#   OPTION: both
#STARparBAM="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM"
#STARparBAM="--outSAMtype BAM Unsorted --quantMode TranscriptomeSAM"
#   OPTION: Gene level counts
#STARparBAM="--outSAMtype BAM Unsorted --quantMode GeneCounts"
#STARparBAM="--outSAMtype BAM Unsorted"
# STAR parameters: strandedness, affects bedGraph (wiggle) files and XS tag in BAM 

case "$dataType" in
str_SE|str_PE)
	#OPTION: stranded data
	STARparStrand=""
	STARparWig="--outWigStrand Stranded"
	;;
#OPTION: unstranded data
unstr_SE|unstr_PE)
	STARparStrand="--outSAMstrandField intronMotif"
	STARparWig="--outWigStrand Unstranded"
	;;
esac

case "$fqType" in
readgroup)
	readGroups=`cat $readGroups`
	STARoutSAMattrRGline="--outSAMattrRGline $readGroups"
	;;
pairs)
	readGroups=`cat $readGroups | head -1 | cut -f 2-`
	STARoutSAMattrRGline="--outSAMattrRGline $readGroups"
	;;
esac


case "$outBamType" in
sorted)
	STARparBAM="--outSAMtype BAM SortedByCoordinate" 
	;;
unsorted)
	STARparBAM="--outSAMtype BAM Unsorted"
	;;
esac
	



# Metadata BAM comments
echo -e '@CO\tREFID:Homo_sapiens_assembly19
@CO\tANNID:gencode.v19' > comments.txt



# STAR parameters: metadata
STARparsMeta="--outSAMheaderCommentFile comments.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate"


###### STAR command
echo "$STAR $STARparCommon $STARparZip $STARparReadLength $STARparRun $STARparBAM $STARparStrand $STARparsMeta $STARoutSAMattrRGline"
$STAR $STARparCommon $STARparZip $STARparReadLength $STARparRun $STARparBAM $STARparStrand $STARparsMeta $STARoutSAMattrRGline
