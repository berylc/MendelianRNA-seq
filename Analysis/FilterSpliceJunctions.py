
import sys ,  argparse , re , itertools
from operator import add
from colorama import init ,  Fore


def filter(args):
	#Make necessary gene lists if specified 
	if args.add_OMIM:
		dOMIM  =  make_omim_dic(args.OMIM_file)
	if args.genes:
		gene_list  =  make_gene_list(args.genes) 

	with open(args.splice_file) as inp:
		for line in inp:
			line  =  line.strip()
			if "Start" in line: 
				continue
			gene , gene_type , junction , n_times , n_samps , all_samps  =  line.split("\t")[0:6]
			n_times ,  n_samps  =  int(n_times) , int(n_samps)
			
			#Exclude junctions observed once in each individual
			if n_times  ==  n_samps: 
				continue
			#Filter junctions by number of samples and read support 
			if n_samps > args.n_samples: 
				continue
			if n_times <=  args.n_read_support:
				continue
			#If gene list specified, filter other genes 
			if args.genes:
				if gene not in gene_list:
					continue 

			#Filter out sample:read pairs if the number of reads in the sample is less than the given value
			filtered_samps  =  []
			kept_samps  =  []
			kept_pairs  =  []
			samp_read_pairs  =  all_samps.split(",")
			if args.filter_n_reads:
				for samp in samp_read_pairs:
					read , sample  =  samp.split(":")
					if int(read) < args.filter_n_reads:
						n_samps  =  n_samps - 1
						n_times  =  n_times - int(read)
						filtered_samps.append(samp)
					else:
						kept_samps.append(sample)
						kept_pairs.append(samp)
			if n_times <= 0 or n_samps  == 0: 
				continue	
			
			new_line  =  [gene ,  gene_type , junction ,  str(n_times) ,  str(n_samps) ,  ",".join(kept_pairs)]
			
			#If a given sample is specified, filter junctions in which the sample does not have the highest read support
			highestsamp  =  kept_pairs[-1]
			if args.sample_with_highest_read_support:
				if args.sample_with_highest_read_support not in highestsamp:
					continue

			#If specified, remove all junctions in where the highest read support is not (x args.times_last) that second highest read support
			if args.times_last:
				try: 
					nexthighest  =  kept_pairs[-2]
					if int(nexthighest.split(":")[0]) * float(args.times_last) > int( highestsamp.split(":")[0])  : 
						continue
				except IndexError:
					pass

			#If specified, remove junctions which are observed highest in GTEx samples
			if args.exclude_gtex_junctions:
				if re.match("^G[0-9]+" , highestsamp.split(":")[1]): continue
		
			#If splice_file has been normalized, apply parameters 
			if args.include_normalized:
				tag , normalized_counts =  line.split("\t")[6:8]
				new_line.insert(4 , tag)	
				if args.normalized_times_last and args.normalized_minimum:
					new_line  =  apply_normalized(new_line , normalized_counts , kept_samps , norm_times_last = args.normalized_times_last , normalized_min = args.normalized_minimum)
				elif args.normalized_times_last:
					new_line  =  apply_normalized(new_line , normalized_counts , kept_samps , norm_times_last = args.normalized_times_last)
				elif args.normalized_minimum:
					new_line  =  apply_normalized(new_line , normalized_counts , kept_samps , normalized_min = args.normalized_minimum)
				else:
					new_line  =  apply_normalized(new_line , normalized_counts , kept_samps)

				if new_line == "-":
					continue

			#If a given sample is specified, filter junctions in which the sample does not have the highest normalized read support
				if args.sample_with_highest_normalized_read_support:
					highest_normalized_samp  =  new_line[0]
					if args.sample_with_highest_normalized_read_support !=  highest_normalized_samp:
						continue
			
			#Add OMIM information 
			if args.add_OMIM:
				OMIMinfo  =  add_OMIM(gene , dOMIM)
				if args.only_OMIM:
					if OMIMinfo  == "-":
						continue
				for col in OMIMinfo: 
					new_line.append(col)

		
			line_to_print  =   "\t".join(new_line)
			if args.print_simple:
				print line_to_print
			else:
				A = re.sub(r'([,|\t][\d|\.]+:[-|\w|\.]+\t)', Fore.RED + r'\1' + Fore.RESET, line_to_print)
				B = re.sub(r'([,|.|\d][\.|\d]+:[\w|\-|A-Z|\s|\.]+)$', Fore.GREEN + r'\1' + Fore.RESET, A)

				if args.add_OMIM:
					B = re.sub(r'([,|.|\d][\.|\d]+:[\w|\-|A-Z|\s|\.]+\t)', Fore.GREEN + r'\1' + Fore.RESET, A)

				print B


def apply_normalized(new_l , norm_counts , k_samps , norm_times_last = None , normalized_min = None , normalizedonlysample = None):
	if norm_counts  ==  "-":
		new_l.insert(0 , "-")
		new_l.append("-")
		return new_l
	normalized_kept_counts = []
	for pair in norm_counts.split(","):
		norm_val , norm_samp  =  pair.split(":")
		if norm_samp in k_samps:
			normalized_kept_counts.append(pair)
	if normalized_kept_counts == []:
		return "-"
	highest_normalized  =  normalized_kept_counts[-1]	
	highest_normalized_samp  =  highest_normalized.split(":")[1]
	highest_normalized_value  =  highest_normalized.split(":")[0]
	new_l.insert(0 , highest_normalized_samp)
	if norm_times_last:
		try:
			second_normalized_value  =  normalized_kept_counts[-2].split(":")[0]
			if "*" in highest_normalized or "*" in second_normalized_value:
				new_l.append(",".join(normalized_kept_counts))		
				second_normalized_value  =  second_normalized_value.split("*")[0]
				highest_normalized_value  =  highest_normalized_value.split("*")[0]
			if float(highest_normalized_value)  >=  float(second_normalized_value)*norm_times_last or len(norm_counts)  == 1:
					new_l.append(",".join(normalized_kept_counts))

			else:
					return "-"
		except IndexError:
			new_l.append(",".join(normalized_kept_counts))
						
	else:
		new_l.append(",".join(normalized_kept_counts))

	if normalized_min:
		if "*" in highest_normalized:
			new_l.append(",".join(normalized_kept_counts))
			return '-'
		if float(highest_normalized_value) < float(normalized_min):
			return '-'
	return new_l


def add_OMIM(gene , omimdic):
	if gene in omimdic.keys():
		info  =  omimdic[gene]
		if len(info)  == 1: 
			return info[0]	
		else: 
			return reshape_omim_info(info)
	else:
		return "-"


def make_omim_dic(omimtable):
	d =  {}
	with open(omimtable) as inp:
		for line in inp:
			genes ,  hgnc_synonyms ,  hgnc_genes ,  phenotype ,  phenotype_inheritance ,  gene_mim_number ,  phenotye_mim_number ,  chrom ,  comments  =  line.strip().split("\t") 
	 		gene_list   =  genes.split("|")
			for gene in gene_list:
				if gene not in d:
					d[gene]   =  []
				d[gene].append([phenotype , phenotype_inheritance])	
	return d


def reshape_omim_info(list_of_omim_input):
	iters  =  [iter(x) for x in list_of_omim_input]
	combined  =  iter(it.next() for it in itertools.cycle(iters))
	combinedv2  =   [','.join(each) for each in itertools.izip(combined , combined)]	
	return combinedv2


def make_gene_list(gene_input):
	try:
		f =  open(gene_input,"r")
		return [x.strip() for x in f.readlines()]
	except IOError:
		genes  =  gene_input.split(",")
		return  genes

if __name__ == "__main__":
	parser  =  argparse.ArgumentParser(description  =  '''Filter splice junction data''')	
	
	input_files_group  =  parser.add_argument_group("Input file parameters")
	input_files_group.add_argument('-splice_file' , help  =  'Input splice file' , action  = 'store' , required = True)
	input_files_group.add_argument('-include_normalized' , help = 'Include if the splice file has been run through NormalizeSpliceJunctionValues.py' , action = 'store_true')

	raw_reads_group  =  parser.add_argument_group("Raw read support-based parameters")
	raw_reads_group.add_argument('-n_read_support' , help = 'Minimum number of reads that support the junction in the dataset' , action = 'store' , type = int ,  default = 0)
	raw_reads_group.add_argument('-filter_n_reads' , help = 'Filter all samples in which a splice junction is supported by less than this number of reads' , action = 'store' , type = int , default = 1)
	raw_reads_group.add_argument('-times_last' , help = 'Filter junctions in which the highest read support is lower than this number times the next highest "' , action = 'store')
	
	normalized_read_group  =  parser.add_argument_group("Normalized read support-based parameters")
	normalized_read_group.add_argument('-normalized_times_last' , help = 'Filter junctions in which the highest normalized read support is lower than this number times the next highest' , action = 'store' , type = float)
	normalized_read_group.add_argument('-normalized_minimum' , help = 'Filter all samples in which a the noramalized read support value is less than this', action = 'store' , type = float)
	
	samples_group  =  parser.add_argument_group("Sample-based parameters")
	samples_group.add_argument('-n_samples' , help = 'Maximum number of samples in which a splice junction is seen' ,  action  =  'store' , type = int , default = 1000)
	samples_group.add_argument('-sample_with_highest_read_support' , help = 'Only output junctions which are seen with highest read support in this sample' , action = 'store')
	samples_group.add_argument("-sample_with_highest_normalized_read_support" , help = 'Only output junctions which are seen with highest normalized read support in this sample' , action = 'store')
	samples_group.add_argument('-exclude_gtex_junctions' , help =  'Exclude splice sites where the highest read supoprt is in a GTEx samples' ,  action  = 'store_true')
	samples_group.add_argument("-exclude_sample" , help = "Exclude junctions where this sample has the highest raw read support" , action = 'store')

	
	gene_lists_group  =  parser.add_argument_group("Gene lists")
	gene_lists_group.add_argument('-add_OMIM' , help = 'Add OMIM information' , action = 'store_true')
	gene_lists_group.add_argument('-OMIM_file' , help = 'Parsed OMIM table, updated versions can be obtained from https://github.com/macarthur-lab/gene_lists' , action = 'store' , default = '/humgen/atgu1/fs03/berylc/MacGit/gene_lists/other_data/omim.use.tsv')
	gene_lists_group.add_argument('-only_OMIM' , help  =  'Only output splicing events seen in OMIM genes' ,  action  =  'store_true')
	gene_lists_group.add_argument('-genes' , help = 'Check splicing in a specific list of genes. Submit a file with genes or comma-seperated string' , action = 'store')
	
	misc  =  parser.add_argument_group("Miscellaneous arguments")
	misc.add_argument("-print_simple" , help = "Print without coloring" , action = 'store_true')
	
	args = parser.parse_args()

	filter(args)	
