import sys


def makeUniqSpliceDict(SpliceFile):
	d = {}
	for iLine in SpliceFile:
		elems = iLine.strip().split()
		try:
			gene, gene_type, sample, chrom, spliceStart, spliceEnd, matchedExon, intronLength = elems
			uniqSplice = "%s:%s:%s:%s-%s"%(gene, gene_type,chrom, spliceStart,spliceEnd)	
		except ValueError:
			gene, gene_type, sample, chrom, spliceStart, spliceEnd, matchedExon, intronLength,detail = elems
			uniqSplice = "%s:%s:%s:%s-%s*%s"%(gene, gene_type,chrom, spliceStart,spliceEnd,intronLength)
		if uniqSplice not in d:
			d[uniqSplice] = {}
		if sample not in d[uniqSplice]:	
			d[uniqSplice][sample]=1

		else:
			d[uniqSplice][sample] += 1
	return d

def printSplices(SpliceFile):
	d = makeUniqSpliceDict(SpliceFile) 
#	print "Gene\tChrom\tStart\tEnd\tNTimesSeen\tNSamplesSeen\tSamples:NSeen"
	for key in d:
		inf = key.split(":")
		Gene = inf[0]
		GeneType = inf[1]
		Chrom = inf[2]
		positions = inf[3].split("-")
		Start = positions[0]
		End = positions[1]
		col3 = []
		nSamplesSeen = len(d[key])
		NTimesSeen = sum(d[key].values())
		for pair in d[key]:
			item = "%s:%s"%(pair,d[key][pair]) 
			col3.append(item)
		SamplesNSeen =  ",".join(col3)
		print "\t".join([str(Gene),str(GeneType),str(Chrom),str(Start),str(End),str(NTimesSeen),str(nSamplesSeen),str(SamplesNSeen)])
	

if __name__=="__main__":
	printSplices(sys.stdin)
