library(ggplot2)
library(argparse)
source('./bin.R')
parser  <-  ArgumentParser(description='Compare patient expression values to GTEX samples to check tissue identity/quality')
parser$add_argument('-gtex_rpkm', action='store', help='GTEx RPKM file,  default in place for muscle check', default='../data/GTEx_RPKM_FB_MSCK_SKN_ADBP_TissuePrefExprGenes.txt')
parser$add_argument('-patient_rpkm', action='store', help='Patient RPKM file i.e. *.genes.rpkm.gct file from RNASeQC output')
parser$add_argument('-outfile', action='store', help='File name prefix to plot PCA')
parser$add_argument('-writePCADat', action='store_true', help='Write text file with PC coordinates', default=FALSE)
args  <-  parser$parse_args()

all_gtex <- read.delim(args$gtex_rpkm, header=T, row.names=1, stringsAsFactors=F)
all_patient <- read.delim(args$patient_rpkm, header=T, skip=2, stringsAsFactors=F)
tissue_list = unique(sapply(strsplit(rownames(all_gtex), split="\\."), function(x) x[1]))

#Subset patient RPKM file to tissue-preferentially expressed genes. 
all_patient <- extractGenesOfInterest(all_patient, makecolnames = F, gene_list=tissue_pref_genes)

#Combine patient and GTEXs values and transform
all_tissues  <-  rbind(all_gtex, all_patient)
all_tissues  <-  log2(all_tissues+1)

#This section performs the PCA,  produces a plot of PCS 1-3 and writes out the PCA text file
thrTissPCA <- prcomp(as.matrix(all_tissues),  retx=TRUE, na.action=na.omit, center=TRUE) 

pdfFile = paste(args$outfile, ".pdf", sep="")
pdf(pdfFile)
plot(thrTissPCA, type="l", main="Variance Explained by PCs")

PCADat <- data.frame(thrTissPCA$x)
PCADat$tissue <- sapply(strsplit(rownames(PCADat), split="\\."),  function(x) x[1])

ggplot(PCADat, aes(x=PC1, y=PC2, col=tissue))+geom_point(size=3)+theme_bw()+theme(plot.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=14), legend.title=element_blank(), legend.text=element_text(size=10))+ggtitle("PC1 vs PC2")

ggplot(PCADat, aes(x=PC1, y=PC3, col=tissue))+geom_point(size=3)+theme_bw()+theme(plot.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=14), legend.title=element_blank(), legend.text=element_text(size=10))+ggtitle("PC1 vs PC3")

ggplot(PCADat, aes(x=PC2, y=PC3, col=tissue))+geom_point(size=3)+theme_bw()+theme(plot.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=14), legend.title=element_blank(), legend.text=element_text(size=10))+ggtitle("PC2 vs PC3")

dev.off()

if (args$writePCADat){
  PCADatFiletoWrite = PCADat[PCADat[, cols+1]=="Patient muscle", ]
  PCADatFiletoWrite = PCADatFiletoWrite[order(rownames(PCADatFiletoWrite)), 1:5]
  Name = paste(args$outfile, ".PCA.txt", sep="")
  write.table(PCADatFiletoWrite,  file = Name,  quote=F,  row.names = T,  col.names=T,  sep="\t") }

#This section performs hiearchichal clustering to output text saying whether your sample clusters closely with muscle. 
#It is just another a way to look at the data,  but the PCA is more useful than this text output. 

#For this you need to specify the number of clusters you expect based on the GTEx data. 
#With skin-lower leg,  skin-subprapubic,  skeletal-muscle,  adipose-subcutaneous and fibroblasts we expect 3.
nClus = 3

hclustfunc  <-  function(x) hclust(x,  method="average")
distfunc  <-  function(x) as.dist(1-cor(t(x), method = "spearman"))
d  <-  distfunc(all_tissues)
fit  <-  hclustfunc(d)
 
allClusters <- cutree(fit, k=nClus)
for (clusnum in 1:nClus){
  clus <-  allClusters[allClusters==clusnum]
  clusnames <- sapply(strsplit(names(clus), split="\\."),  function(x) x[1])
  clusnames = unique(clusnames)
  if ("muscle" %in% clusnames){
    cat(paste("\n", "Patient samples that cluster with GTEx muscle samples"))
    for (elem in sort(clusnames)){
      if (elem=="muscle"){next}
      if (substr(elem,  1,  1)=="X"){cat(paste("\n", substr(elem, 2, nchar(elem))))
    } 
      else{cat(paste("\n", elem))
    }
    }
    }
   
  if (!"muscle" %in% clusnames && length(setdiff(clusnames, tissue_list))!=0){
    cat(paste("\n\n", "Patient samples that do not cluster tightly with muscle:"))
     
    for (eachSamp in setdiff(clusnames, tissue_list)){cat(paste("\n", eachSamp))                                                   
    cat(paste("\n", "PCA coordinates are:", paste(PCADat[eachSamp, c("PC1", "PC2")], collapse=", "), '\n'), "\n") 
                                                     
     }
     }
   }
   
