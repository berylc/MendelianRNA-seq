library(ggplot2)
library(argparse)
source('./bin.R')
parser  <-  ArgumentParser(description='Compare patient expression values to GTEX samples to check tissue identity/quality')
parser$add_argument('-gtex_rpkm', action='store', help='GTEx RPKM file,  default in place for muscle check', default='../data/gtex_expression_tissue_preferential_fibs.msck.skn.adbp.txt')
parser$add_argument('-patient_rpkm', action='store', help='Patient RPKM file i.e. *.genes.rpkm.gct file from RNASeQC output')
parser$add_argument('-out_file', action='store', help='File name prefix to plot PCA')
parser$add_argument('-writePCADat', action='store_true', help='Write text file with PC coordinates', default=FALSE)
args  <-  parser$parse_args()

all_gtex <- read.delim(args$gtex_rpkm, header=T, row.names=1, stringsAsFactors=F)
all_patient <- read.delim(args$patient_rpkm, header=T, skip=2, stringsAsFactors=F)
names(all_patient) <- gsub("\\.","-",names(all_patient))

tissue_list = unique(gsub("\\.[0-9]*$","",rownames(all_gtex)))

#Subset patient RPKM file to tissue-preferentially expressed genes. 
all_patient <- extractGenesOfInterest(all_patient, gene_list=tissue_pref_genes)
 
#Combine patient and GTEXs values perform PCA
genes_in_both <- intersect(names(all_patient), names(all_gtex))
all_tissues <- rbind(subset(all_gtex, select = genes_in_both), subset(all_patient, select = genes_in_both))

PCADat <- performPCA(all_tissues, "tissue",expected_annotations =tissue_list )

#Produces a plot of PCS 1-3
plot_theme <- theme_bw()+theme(plot.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=14), legend.title=element_blank(), legend.text=element_text(size=10),panel.grid.major=element_blank())
colours <- scale_colour_manual(values=c("#FF8C69","#EEB422", "#009900","#33CC33","#1E90FF","#104E8B"))

pdfFile = paste(args$out_file, ".pdf", sep="")
pdf(pdfFile,width=10 , height = 6)
plot(PCA, type="l", main="Variance Explained by PCs")

ggplot( PCADat, aes( x = PC1, y = PC2, col = tissue ))+geom_point( size=3 ) +
              plot_theme + colours +
              xlab( paste( "PC1 (", signif( PoV[1]*100, 2 ), "%)", sep = "" ))+ylab( paste (" PC2 (",signif( PoV[2]*100, 2 ), "%)", sep = ""))

ggplot(PCADat, aes(x=PC1, y=PC3, col=tissue))+geom_point(size=3) +
              plot_theme + colours +
              xlab( paste( "PC1 (", signif( PoV[1]*100, 2 ), "%)", sep = "" ))+ylab( paste (" PC3 (",signif( PoV[3]*100, 2 ), "%)", sep = ""))


ggplot(PCADat, aes(x=PC2, y=PC3, col=tissue))+geom_point(size=3)+
              plot_theme + colours +
              xlab( paste( "PC2 (", signif( PoV[2]*100, 2 ), "%)", sep = "" ))+ylab( paste (" PC3 (",signif( PoV[3]*100, 2 ), "%)", sep = ""))


dev.off()

#Write out the PCA information to a text file
if (args$writePCADat){
  PCADatFiletoWrite = PCADat[PCADat[, ncol(PCADat)]=="Patient", ]
  PCADatFiletoWrite = PCADatFiletoWrite[order(rownames(PCADatFiletoWrite)), 1:5]
  Name = paste(args$out_file, ".PCA.txt", sep="")
  write.table(PCADatFiletoWrite,  file = Name,  quote=F,  row.names = T,  col.names=T,  sep="\t") }

#The following section performs hiearchichal clustering to output text indicated whether the input sample clusters closely with muscle. 
#It is just another a way to look at the data

#For this you need to specify the number of clusters you expect based on the GTEx data. 
#With skin-lower leg,  skin-subprapubic,  skeletal-muscle,  adipose-subcutaneous and fibroblasts we expect 3 (you can see this from the PCA too)
nClus = 3

hclustfunc  <-  function(x) hclust(x,  method="average")
distfunc  <-  function(x) as.dist(1-cor(t(x), method = "spearman"))
d  <-  distfunc(all_tissues)
fit  <-  hclustfunc(d)
 
allClusters <- cutree(fit, k=nClus)
for (clusnum in 1:nClus){
  clus <-  allClusters[allClusters==clusnum]
  clusnames <- unique(gsub("\\.[0-9]*$","",names(clus)))
  if ("Muscle.Skeletal" %in% clusnames){
    cat(paste("\n", "Patient samples that cluster with GTEx muscle samples"))
    for (elem in sort(clusnames)){
      if (elem=="Muscle.Skeletal"){next}
      cat(paste("\n", elem))
    }
  }
  
  if (!"Muscle.Skeletal" %in% clusnames && length(setdiff(clusnames, tissue_list))!=0){
    cat(paste("\n\n", "Patient samples that do not cluster tightly with muscle:"))
    
    for (eachSamp in setdiff(clusnames, tissue_list)){cat(paste("\n", eachSamp))                                                   
      cat(paste("\n", "PCA coordinates are:", paste(PCADat[eachSamp, c("PC1", "PC2")], collapse=", "), '\n'), "\n") 
      
    }
  }
}

