#!/bin/python 
import time
import argparse
import os
import subprocess
import shutil 
import sys
def main(args):
	#Validate files and extract relevant information 
	if args.bam.endswith("bam"):
		if args.data_type=="RNA":
			sampname = args.bam.strip().split("/")[-1].strip(".sorted.deduped.rg_added.split.bam")
		else: 
			sampname = args.bam.strip().split("/")[-1].strip(".bam")
	else:
		raise Exception("Bam file must end in .bam")

	#Define reference files and region to extract reads from	
	if args.gene =='TTN':
		trip_region = "2:179517931-179528342"
		if args.region   == "1": reference = args.ttn_trip1_reference
		elif args.region == "2": reference = args.ttn_trip2_reference
		elif args.region == "3": reference = args.ttn_trip3_reference
	if args.gene == 'NEB':
		trip_region  = "2:152431602-152470936"
		if args.region   == "1"  : reference = args.neb_trip1_reference
		elif args.region == "2": reference = args.neb_trip2_reference
	refname = reference.replace(".fasta","")
	
	#Make directory 
	new_dir = "%s/%s"%(os.getcwd(),sampname)
	if os.path.exists(new_dir):
		print "There seems to be a directory for %s, removing it now to make a new one\n"%sampname
		shutil.rmtree(new_dir)
	
	os.makedirs(new_dir)

	#Run 
	print "1) Extracting reads from %s region" %(args.gene)
	basename = "%s/%s_%s_%s_%s_Region."%(new_dir,sampname,args.data_type,args.gene,args.region)
	samtools_query = ["samtools","view",args.bam,trip_region,"-bh","-o","".join([basename,"bam"])]
	print " ".join(samtools_query),"\n"
	subprocess.check_output(samtools_query)

	print "2) Bam To Fastq"
	bedtools_query = ["bedtools","bamtofastq","-i","".join([basename,"bam"]),"-fq","".join([basename,"fastq"])]
	print " ".join(bedtools_query),"\n"
	subprocess.check_output(bedtools_query)

	print "3) Aligning with BWA. All samples will have the same toy read group information.\n"
	RG =args.read_group.replace("SAMPLENAME",sampname)
	RG = RG.replace("\t","\\t")
	bwa_query = ["bwa","mem","-R",RG,refname,"".join([basename,"fastq"])]
	sam_file = open("".join([basename,"sam"]),"w")
	subprocess.call(bwa_query,stdout=sam_file)
	sam_file.close()
	print_line_of_sam("".join([basename,"sam"]))

	print "\n4) Sorting SAM and indexing\n"
	sort_sam_query = ["java","-jar","-Xmx16g",args.picard_jar,"SortSam","INPUT=%s"%("".join([basename,"sam"])),"OUTPUT=%s"%("".join([basename,"sorted.bam"])),"SORT_ORDER=coordinate","VALIDATION_STRINGENCY=LENIENT"]
	subprocess.check_output(sort_sam_query)

	index_query = ["samtools","index","".join([basename,"sorted.bam"])]
	subprocess.check_output(index_query)

	print "\n 5) Performing variant calling with Haplotype Caller for %s"%(args.data_type)
	if args.data_type =="WES" or args.data_type=="WGS":
		stand_call_conf = "30.0" 
		stand_emit_conf = "10.0"
	if args.data_type == "RNA":
		stand_call_conf = "20.0" 
		stand_emit_conf = "20.0"

	haplotype_caller_query = ["java","-jar","-Xmx4g",args.gatk_jar,"-T","HaplotypeCaller","-R",reference,"-I","".join([basename,"sorted.bam"]),"-o","".join([basename,"vcf.gz"]),"-ploidy","6","-stand_call_conf",stand_call_conf,"-stand_emit_conf",stand_emit_conf ]
	subprocess.check_output(haplotype_caller_query)

	print "\n 6) Running VEP annotation on the mini vcf"
	vep_query = ["/broad/software/free/Linux/redhat_5_x86_64/pkgs/perl_5.10.1/bin/perl","/humgen/atgu1/fs03/konradk/vep/ensembl-tools-release-80/scripts/variant_effect_predictor/variant_effect_predictor.pl","--everything","--vcf","--allele_number"," --no_stats","--cache","--offline","--dir", "/humgen/atgu1/fs03/weisburd/xbrowse/data/vep_cache/","--force_overwrite","--cache_version","80","--fasta","/humgen/atgu1/fs03/birnbaum/loftee-dev/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa","--assembly","GRCh37","--fork","4","-i","".join([basename,"vcf.gz"]),"-o","".join([basename,"vep.vcf.gz"])]
	subprocess.check_output(vep_query)

	if not args.keep_all_files:
		print "\n Deleting the fastq,sam and bam files to save space"
		os.remove("".join([basename,"sam"]))
		os.remove("".join([basename,"fastq"]))
		os.remove("".join([basename,"bam"]))
	


def print_line_of_sam(sam):
	with open(sam) as inp:
		for line in inp:
			if line.strip().startswith("@"):
				continue
			else:
				print "\nFirst line of realigned sam file is"
				print line
				break




class ForceIOStream:
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()
        if not self.stream.isatty():
            os.fsync(self.stream.fileno())

    def __getattr__(self, attr):
        return getattr(self.stream, attr)







if __name__== "__main__":
	parser = argparse.ArgumentParser(description = '''Remap the triplicated region of NEB or TTN''')
	query = parser.add_argument_group('Input options')
	query.add_argument('-bam',help='Bam to realign',required=True)

	query.add_argument('-data_type',help='RNA,WES or WGS',choices=['WES','RNA','WGS'],required=True)
	query.add_argument('-gene',help='Gene to realign the triplicate region',choices=['TTN','NEB'],required=True)
	query.add_argument('-region',help='Triplicate region to realign to',choices=["1","2","3"],required=True)
	references = parser.add_argument_group('Reference options')
	references.add_argument('-neb_trip1_reference',default='/humgen/atgu1/fs03/berylc/MuscDisease/tmp/triplicateRegion/v2/data/references/trip1ref/hg19.onlytrip1.fasta')
	references.add_argument('-neb_trip2_reference',default='/humgen/atgu1/fs03/berylc/MuscDisease/tmp/triplicateRegion/v2/data/references/trip2ref/hg19.onlytrip2.fasta')
	references.add_argument('-ttn_trip1_reference',default='/humgen/atgu1/fs03/berylc/MuscDisease/tmp/triplicateRegion/v3TTN/reference/TTNtrip1ref/hg19.ttn.onlytrip1.fasta')
	references.add_argument('-ttn_trip2_reference',default='/humgen/atgu1/fs03/berylc/MuscDisease/tmp/triplicateRegion/v3TTN/reference/TTNtrip2ref/hg19.ttn.onlytrip2.fasta')
	references.add_argument('-ttn_trip3_reference',default='/humgen/atgu1/fs03/berylc/MuscDisease/tmp/triplicateRegion/v3TTN/reference/TTNtrip3ref/hg19.ttn.onlytrip3.fasta')
	misc = parser.add_argument_group('Miscellaneous arguments')
	misc.add_argument('-read_group',default= "@RG\tID:C412H.5\tLB:Solexa-316447\tSM:SAMPLENAME\tPL:illumina",help="Read group information to add to realigned bam")
	misc.add_argument('-picard_jar',default="/seq/software/picard/current/bin/picard.jar",help="Path to Picard .jar")
	misc.add_argument('-gatk_jar',default="/humgen/gsa-hpprojects/GATK/bin/current/GenomeAnalysisTK.jar",help="Path to GATK jar")
	misc.add_argument('-keep_all_files',action='store_true',help='If specified won\'t delete any output files')
	args = parser.parse_args()
	sys.stdout = ForceIOStream(sys.stdout)
	sys.stderr = ForceIOStream(sys.stderr)
	print "Started on", time.strftime("%d/%m/%Y"), "at", time.strftime("%H:%M")
	main(args)
	print "Done! Finished on", time.strftime("%d/%m/%Y"), "at", time.strftime("%H:%M")
