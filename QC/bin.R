tissue_pref_genes <- read.delim("../data/tissue_preferential_genes_fibs.msck.skn.adbp.txt",header=T,stringsAsFactors=F)$ensg_id
sex_biased_genes <- read.delim("../data/sex_biased_genes.txt",header=T,stringsAsFactors = F,strip.white=T)


extractGenesOfInterest<-function(rpkm_file,gene_list,patient_col=3){
    gene_list <- sapply(strsplit(gene_list, split="\\."), function(x) x[1])
    rpkm_file$Name <- sapply(strsplit(rpkm_file$Name, split="\\."), function(x) x[1])
    
    rpkm_file<-subset(rpkm_file, Name %in% gene_list, select = -c(Description))
    
    rownames(rpkm_file)<-rpkm_file$Name
    rpkm_file <- subset(rpkm_file, select = -c(Name))
    rpkm_file <- data.frame(t(rpkm_file),stringsAsFactors=F)
    return(rpkm_file)
}


performPCA <- function(rpkm_file, annotation_column, expected_annotations){
  rpkm_file  <-  log2(rpkm_file+1)
  PCA <<- prcomp(as.matrix(rpkm_file),  retx=TRUE, na.action=na.omit, center=TRUE)

  PoV <<- PCA$sdev^2/sum(PCA$sdev^2)
  PCADat <- data.frame(PCA$x)

  PCADat[,annotation_column] <- gsub("\\.[0-9]*$","",rownames(PCADat))

  PCADat[!PCADat[,annotation_column] %in% expected_annotations, annotation_column] <- "Patient"

  return(PCADat)
}

getXISTAvgYChromExpr <- function(rpkm_file, chrom_y_genes, annotation_column){
  chrom_y_genes <- sapply(strsplit(chrom_y_genes, split="\\."), function(x) x[1])
  chrom_y_genes <- intersect(names(rpkm_file), chrom_y_genes)
  
  XIST <- subset(rpkm_file, select = ENSG00000229807 )
  chrom_y_genes <- subset(rpkm_file, select = chrom_y_genes )
  chrom_y_genes$avg_y_chrom <- rowMeans(chrom_y_genes)
  
  xist_y_chrom <- merge(XIST, subset(chrom_y_genes, select= avg_y_chrom), by='row.names')
  
  xist_y_chrom[,"sex"] <- gsub("\\.[0-9]*$","",xist_y_chrom$Row.names)
  xist_y_chrom[!xist_y_chrom$sex %in% c("male","female"), "sex"] <- "Patient"
  
  names(xist_y_chrom) <- c("sample","XIST","avg_y_chrom","sex")
  return(xist_y_chrom)
  
}