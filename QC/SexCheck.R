library(argparse)
library(ggplot2)
library(plyr)
library(dplyr)
source('./bin.R')
parser <- ArgumentParser(description='Compare patient expression values to GTEX samples to confirm sex. Currently set up for muscle tissue.')
parser$add_argument('-patient_rpkm',help='Patient RPKM file (i.e. *.genes.rpkm.gct file from RNASeQC output)')
parser$add_argument('-gtex_sex_biased',help='Precomputed sex biased expression file from GTEX',default='../data/gtex_expression_sex_biased.txt')
parser$add_argument('-out_file',help='File name prefix for output file')
args <- parser$parse_args()

gtex_sex_biased_rpkm<-read.delim(args$gtex_sex_biased,stringsAsFactors=F, row.names=1)
patients_rpkm<-read.delim(args$patient_rpkm,header=T, stringsAsFactors = F,skip=2)

patients_rpkm<- extractGenesOfInterest(patients_rpkm, sex_biased_genes$gene_id)

genes_in_both <- intersect(names(gtex_sex_biased_rpkm), names(patients_rpkm))
all_sex_biased <- rbind(subset(gtex_sex_biased_rpkm, select = genes_in_both), subset(patients_rpkm, select = genes_in_both))

PCADat <- performPCA(all_sex_biased, "sex",expected_annotations =c("male","female") )

chrom_y_genes <-  subset(sex_biased_genes, chr == "chrY")$gene_id
all_xist_y <- getXISTAvgYChromExpr(all_sex_biased,chrom_y_genes )

plot_theme <- theme_bw()+theme(plot.title=element_text(size=20), axis.text=element_text(size=15), axis.title=element_text(size=14), legend.title=element_blank(), legend.text=element_text(size=10),panel.grid.major=element_blank())
colours <- scale_colour_manual(values=c("#700D4F","#EAA8D4","#87B2DD"))


pdfFile = paste(args$out_file, ".pdf", sep="")
pdf(pdfFile, width =8, height =6)

plot(PCA,type="l",main="Variance Explained by PCs")

ggplot( PCADat %>% arrange(sex), aes( x = PC1, y = PC2, col = sex))+geom_point( size=3  ) + 
  plot_theme + colours +
  xlab( paste( "PC1 (", signif( PoV[1]*100, 2 ), "%)", sep = "" ))+ylab( paste (" PC2 (",signif( PoV[2]*100, 2 ), "%)", sep = ""))

ggplot( all_xist_y %>% arrange(sex), aes( x = XIST, y = avg_y_chrom, col = as.factor(sex)))+geom_point(size = 3)+
  plot_theme + colours + 
  xlab("XIST expression") + ylab("Average of expression of sex-biased genes on y chromosome")

dev.off()

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) as.dist(1-cor(t(x),method = "spearman"))
d <- distfunc(all_sex_biased)
fit <- hclustfunc(d)
allClusters<-cutree(fit, k = 2)


for (clusnum in 1:2){
  clus <-  allClusters[allClusters==clusnum]
  cluster <- unique(gsub("\\.[0-9]*$","",names(clus)))
  cat('\n')
  if('male' %in% cluster){
    cat("\nMale patients are:")
    for(elem in names(clus)){
      if(grepl("male",elem)){next}
      cat(paste("\n", elem))
      }
  }
  if('female' %in% cluster){
    cat("\nFemale patients are:")
    for(elem in names(clus)){
      if(grepl("female",elem)){next}
      cat(paste("\n", elem))
      }
  }
  }

    