#!/bin/bash/py
import sys
import argparse
import gzip
import re
import itertools
from operator import add

def make_annotated_junction_set(f,chrom_col,start_col,stop_col):
	'''From an annotations file, extracts all all annotated junctions with 1 bp flank on either side. 
		These are junctions as in 1:1221-1345. So these are annotated splice junctions'''
	s=set()
	with f  as inp:
		for junction in inp:
			fields = junction.strip().split("\t") 
			chrom = fields[chrom_col]
			start= fields[start_col]
			stop = fields[stop_col]
			chrom = chrom.strip("chr")
			start_flanks = get_flank(start)
			stop_flanks = get_flank(stop)
			st = tuple(("%s:%s-%s"%(chrom,start_flanks[ind],stop_flanks[ind])) for ind in range(0,3))
			stadd1 = "%s:%s-%s"%(chrom,start_flanks[0],stop_flanks[2])				#Need to think about these flanks
			stadd2 = "%s:%s-%s"%(chrom,start_flanks[2],stop_flanks[0])
			s.update(st,stadd1,stadd2)
	return s


def get_flank(pos):
	pos = int(pos)-1
	return [pos-1,pos,pos+1]


def sort_floats(pair_list,type="float"):
	'''Sorting SampleName:Nread by read count to place the highest at the end and colour t'''
 	count_dict = {}
 	sorted_floats = []
	for pair in pair_list:
		sample,n_reads = pair.split(":")
		if "*" in n_reads: #If the annotated junction doesn't align reads currently not able to detect them efficiently 
			sorted_floats.insert(0,"%s:%s"%(n_reads,sample))
			continue
		if type=="int":
			n_reads = int(n_reads)
		else:
			n_reads = float(n_reads)
		if n_reads not in count_dict:
			count_dict[n_reads] = [sample]
		else:
			count_dict[n_reads].append(sample)
	for counts in sorted(count_dict.keys()):
		for elem in count_dict[counts]:
				sorted_floats.append("%s:%s"%(counts,elem))
	return sorted_floats

def get_unannotated_junctions(splicefile,annotated,chrom_col,start_col,stop_col):
	with open(splicefile) as inp:
		for spliceline in inp:			
			fields = spliceline.strip().split("\t")
			chrom = fields[2]
			start= fields[3]
			stop = fields[4]
			chrom = chrom.strip("chr")
			#gene,chrom,start,stop,ntimes,nsamp,samptimes = spliceline.strip().split("\t")
			q = "%s:%s-%s"%(chrom,start,stop)
			if q in annotated: 
				continue
			else: 
				print spliceline.strip()


def get_annotated_counts(splicefile,annotated):
	''' '''
	splice_dict = {}
	with open(splicefile) as inp:  
		for spliceline in inp: 
			if "Chrom" in spliceline : 
				continue
			gene,gene_type,chrom,start,stop,ntimes,nsamp,samptimes = spliceline.strip().split("\t")[0:8]
			q = "%s:%s-%s"%(chrom,start,stop)
			#If splice junction is not annotated
			if q not in annotated: 
				continue
			'''If junction is annotated, note the number of reads that align to that exon-intron junction in all samples that carry it
			If your splice junction is annotated and the junction is ~ 2:100-250 Besse:10, Beryl 20 then the dictionary you create is:
			{'2:200':{'Beryl':10,'Besse':10}  ,  2:100: {'Beryl':10,'Besse':10}}'''
			for pos in [start,stop]:
				pos = "%s:%s"%(chrom,pos) 
				if pos not in splice_dict:
					splice_dict[pos] = {} #Initialize 
					for pair in samptimes.split(","):
						pair_sample = pair.split(":")[0]
						pair_ntimes = pair.split(":")[1]
						if pair_sample not in splice_dict[pos]:
							splice_dict[pos][pair_sample]=int(pair_ntimes)
				if pos in splice_dict: 												#If the position is in dictionary, then update on a per sample basis to get the max splice junction reads
					current_counts = splice_dict[pos]
					for pair in samptimes.split(","):
						pair_sample = pair.split(":")[0]
						pair_ntimes = pair.split(":")[1]
						current_count_sample = current_counts.get(pair_sample,0)
						if int(pair_ntimes) > int(current_count_sample):
							current_counts[pair_sample] = int(pair_ntimes)
	return splice_dict



def normalize_counts(splicefile,annotated_counts): 
	with open(splicefile) as inp:
	 	for spliceline in inp:
	 		if "Chrom" in spliceline:
	 			continue
			normalizedCol = []
			gene,gene_type,chrom,start,stop,ntimes,nsamp,samptimes = spliceline.strip().split("\t")  
			#For printing later 
			full_junction = "%s:%s-%s"%(chrom,start,stop)
			samptimes_sorted = sort_floats(samptimes.split(","),type="int")
			base_start = "%s:%s"%(chrom,start) 
			base_stop= "%s:%s"%(chrom,stop)
			annotated_start = annotated_counts.get(base_start)
			annotated_stop = annotated_counts.get(base_stop)
			if annotated_start and  annotated_stop:													#Canonical splicing and exon skipping
				tag = "Both annotated"
				for pair  in samptimes.split(","):
					pair_sample = pair.split(":")[0]
					pair_ntimes = int(pair.split(":")[1])
					annotated_start_sample = annotated_start.get(pair_sample,0)
					annotated_stop_sample = annotated_stop.get(pair_sample,0)
					denominator  = max(annotated_start_sample,annotated_stop_sample)
					try:
						normalized = round(float(pair_ntimes)/int(denominator),3)
					except ZeroDivisionError:
						normalized = "".join([str(pair_ntimes),"*"])
					normalizedCol.append("%s:%s"%(pair_sample,normalized))
				normalizedCol_sorted = sort_floats(normalizedCol)
				line_to_print = "\t".join([gene,gene_type,full_junction,ntimes,nsamp,",".join(samptimes_sorted) ,tag,",".join(normalizedCol_sorted)])
				print line_to_print
				continue
			if annotated_start or annotated_stop:										#When one end of splice junction is annotated: exon extension and intron inclusion
				tag = "One annotated"
				annotated_one = [x for x in [annotated_start,annotated_stop] if x is not None][0]
				for pair  in samptimes.split(","):
					pair_sample = pair.split(":")[0]
					pair_ntimes = int(pair.split(":")[1])
					denominator = annotated_one.get(pair_sample,0) #If there are no reads at the annotated splice junction tag it
					try:
						normalized = round(float(pair_ntimes) / int(denominator),3)
					except ZeroDivisionError:
						normalized = "".join([str(pair_ntimes),"*"])
					normalizedCol.append("%s:%s"%(pair_sample,normalized))
				normalizedCol_sorted = sort_floats(normalizedCol)
				line_to_print = "\t".join([gene,gene_type,full_junction,ntimes,nsamp,",".join(samptimes_sorted) ,tag,",".join(normalizedCol_sorted)])
				print line_to_print
				continue
			if annotated_stop is None and annotated_stop is None:
				tag = "Neither annotated"
			 	line_to_print = "\t".join([gene,gene_type,full_junction,ntimes,nsamp,",".join(samptimes_sorted),tag,"-"])
				print line_to_print

			

def main(args):
	if args.gzipped: 
		f = gzip.open(args.transcript_model)
	else: 
		f = open(args.transcript_model)
	s = make_annotated_junction_set(f,args.chrom_col,args.start_col,args.stop_col)
	if args.normalize:
		annotated_counts = get_annotated_counts(args.splice_file,s)
		normalize_counts(args.splice_file,annotated_counts)
	if args.getunannotated:
		get_unannotated_junctions(args.splice_file,s,args.chrom_col,args.start_col,args.stop_col)







if __name__=="__main__":
	parser = argparse.ArgumentParser(description = '''Get unannotated junctions from splice file''')
	parser.add_argument('-transcript_model',help="Transcript model of canonical splicing, e.g. gencode v19. Default is /humgen/atgu1/fs03/berylc/MuscDisease/bin/references/gencode.comprehensive.splice.junctions.txt",action='store',default = "/humgen/atgu1/fs03/berylc/MuscDisease/bin/references/gencode.comprehensive.splice.junctions.txt")
	parser.add_argument('-splice_file',help='Splice junction file to filter')
	parser.add_argument('-gzipped', help='Add if sjout file is gzipped',action='store_true')
	parser.add_argument('-chrom_col',help='Chromosome column',type=int,default=0)
	parser.add_argument('-start_col',help='Junction start column',type=int,default=1)
	parser.add_argument('-stop_col',help='Junction stop column',type=int,default=2)
	mode_arguments = parser.add_mutually_exclusive_group(required=True)
	mode_arguments.add_argument('--getunannotated',action='store_true',help='From splice junction file, will only return ones not in sjout')
	mode_arguments.add_argument('--normalize',action='store_true',help='Local normalization on splice junction file based on annotated junctions')
	args=parser.parse_args()	
	main(args)
