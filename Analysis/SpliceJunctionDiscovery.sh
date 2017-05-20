#This version (v6) has the sed command to account for softclipped reads

transcriptFile=$1
bamList=$2
baseDir=`pwd`
IFS=$'\n'
groupname=`basename $transcriptFile | sed 's/gencode.v19.genes.filtered.group//'`
echo "Started on" $(date)
echo "Working in directory $baseDir"
echo "Transcript file is $transcriptFile"
echo "Identifying splice junction is $bamList"

echo -e "Gene\tType\tChrom\tStart\tEnd\tNTimesSeen\tNSamplesSeen\tSamples:NSeen" > All.${groupname}.splicing.txt
for line in `cat $transcriptFile`
do
start=`echo $line | cut -f5`
stop=`echo $line | cut -f6`
chrom=`echo $line | cut -f4`
gene=`echo $line | cut -f1`
gene_type=`echo $line | cut -f7`
base=$baseDir/$gene
pos=$chrom":"$start"-"$stop
echo 'processing' $gene
for i in `cat $bamList`
do

sample=`basename $i | sed 's/.sorted.deduped.bam//'`
samtools view ${i} ${pos} | awk '($6 ~ /N/)' |awk '$5==60'| awk 'int($2)< 256' | awk -v sta=$start -v sto=$stop '$4>sta&&$4<sto {print $4,$6}' | cut -d "N" -f 1 | tr 'M' ' ' |  sed -r 's/[0-9]+S//' | awk -v s=$sample -v ge=$gene -v t=$gene_type -v chr=$chrom '{print ge,t,s,chr,$1+$2-1,$1+$2+$3,$2,$3,$4}' >> ${base}.splicing.txt
done
echo $gene 'is complete'
python ./SpliceJunctionSummary.py <  ${base}.splicing.txt >> All.${groupname}.splicing.txt
rm ${base}.splicing.txt
done

echo "Finished on" $(date)


