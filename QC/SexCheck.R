library(argparse)
library(ggplot2)
source("/humgen/atgu1/fs03/berylc/MuscDisease/QC/SexCheck/scripts/sourceSex.R")
parser <- ArgumentParser(description='Sex Check QC')

parser$add_argument('-gtexphenotype',help='GTEx phenotype file, default stored',default='/humgen/atgu1/fs03/berylc/MuscDisease/QC/SexCheck/data/GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt')

parser$add_argument('-patientrpkm',help='Gene-level RPKM file output ')

parser$add_argument('-sexbiasedgenes',help='Sex Biased Genes to differentiate samples on, Mele genelist stored',default='/humgen/atgu1/fs03/berylc/MuscDisease/QC/SexCheck/data/MeleSexBiasedGenelist.txt')

parser$add_argument('-gtexXIST_ychrom',help='Precomputed XIST and Ychrom expression file from GTEX',default='/humgen/atgu1/fs03/berylc/MuscDisease/QC/SexCheck/data/GTExGender_Xist_YchromAvg_expression.txt')

parser$add_argument('-gtex_allSexBiased',help='Precomputed sex biased expression file from GTEX',default='/humgen/atgu1/fs03/berylc/MuscDisease/QC/SexCheck/data/GTExGender_AllSexBiased_expression.txt')


parser$add_argument('-out',help='PDF to write plots to')
args <- parser$parse_args()

sexBiasedGenes<-read.delim(args$sexbiasedgenes,header=T,stringsAsFactors = F,strip.white=T)
sexBiasedGenes_yChrom<-subset(sexBiasedGenes,chr=="chrY")
gtex_phenotypes<-read.delim(args$gtexphenotype,header = T,stringsAsFactors = F)
patients_rpkm<-read.delim(args$patientrpkm,header=T, stringsAsFactors = F,skip=2)

gtex_XIST_yhcrom<-read.delim(args$gtexXIST_ychrom,stringsAsFactors=F)

gtex_allSexBiased<-read.delim(args$gtex_allSexBiased,stringsAsFactors=F)

getSexBiased_xist_y(patients_rpkm,patientcolbegin=3,out="patients",sexBiasedGenes_yChrom)

getcombinedPCAData(gtex_allSexBiased,patients_rpkm,out="all",sexBiasedGenes=sexBiasedGenes)

PCADat_pcs<-data.frame(PCADat_all$x)
addPhenotype(All=PCADat_pcs,"bothPCS",gtex_phenotypes=gtex_phenotypes)
names(WithSex_bothPCS)[ncol(WithSex_bothPCS)]<-"gender"
XIST_Y_patients$gender="Patient"
pdf(args$out)


ggplot(gtex_XIST_yhcrom,aes(x=XIST,y=avgYchrom,col=gender))+geom_point(size=4)+theme_classic()+theme(plot.title=element_text(size=25),axis.text=element_text(size=15),axis.title=element_text(size=20),legend.title=element_blank(),legend.text=element_text(size=20))+ylab("mean expression \n from Y chrom sex biased genes") + xlab("XIST expression")+geom_point(data=XIST_Y_patients,aes(x=XIST,y=avgYchrom),size=4)



plot(PCADat_all,type="l",main="Variance Explained by PCs- GTEX")
ggplot(WithSex_bothPCS,aes(x=PC1,y=PC2,col=gender))+geom_point(size=3)+theme_bw()+theme(plot.title=element_text(size=25),axis.text=element_text(size=15),axis.title=element_text(size=20),legend.title=element_blank(),legend.text=element_text(size=20))+scale_colour_manual(values=c("pink","blue","black"))+ggtitle("PC1 vs PC2")

 WithSex_bothPCS[grep("GTEX",rownames(WithSex_bothPCS)) && WithSex_bothPCS$gender=="Male","IndicatedSex"]<-"GTEx Male"
# WithSex_bothPCS[grep("GTEX",rownames(WithSex_bothPCS)) && WithSex_bothPCS$gender=="Female","IndicatedSex"]<-"GTEx Female"
# WithSex_bothPCS[rownames(WithSex_bothPCS) %in% c("BON_UC219.1_1", "BON_B12.74.1_1","BON_B16.22_1", "BON_UC473_1","BON_B16.19_1", "BON_B14.20_1"),"IndicatedSex"]<-"Patient Male"
# 
# WithSex_bothPCS[rownames(WithSex_bothPCS) %in% c("BON_B14.75.1_1" ,"BON_B12.33.2_1","BON_B09.27.1_1","BON_B14.71.2_1"),"IndicatedSex"]<-"Patient Female"


#WithSex_bothPCS[rownames(WithSex_bothPCS$)]

dev.off()



#Not normalizing at all 
both<-as.matrix(AllSexBiased_RPKM_all)
#x <- scale(both)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) as.dist(1-cor(t(x),method = "spearman"))
d <- distfunc(both)
fit <- hclustfunc(d)
x<-cutree(fit, k = 2)
Clus1<-names(x[x==1])
Clus2<-names(x[x==2])


#Clus1 is all females
determineSex <- function (cluster){
  if(cluster==1){cat("Male patients are:")}
  if(cluster==2){cat("Female patients are:")}
}


Clus1Short<-c()
for(i in Clus1){
  id = strsplit(as.character(i),split="\\.")[[1]]
  subID = paste(id[1],"-",id[2],sep="")
  Clus1Short<-c(Clus1Short,subID)
}
Clus2Short<-c()

for(i in Clus2){
  id = strsplit(as.character(i),split="\\.")[[1]]
  subID = paste(id[1],"-",id[2],sep="")
  Clus2Short<-c(Clus2Short,subID)
}




patientClus2<-Clus2[-grep("GTEX",Clus2)]
patientClus1<-Clus1[-grep("GTEX",Clus1)]

cat("\n")
determineSex(unique(gtex_phenotypes[gtex_phenotypes$SUBJID%in%Clus1Short,"GENDER" ]))
for(elem in patientClus1){
  if(substr(elem, 1, 1)=="X"){cat(paste("\n",substr(elem,2,nchar(elem))))} 
  else{cat(paste("\n",elem))}}
  
cat("\n")
cat("\n")
determineSex(unique(gtex_phenotypes[gtex_phenotypes$SUBJID%in%Clus2Short,"GENDER" ]))

for(elem in patientClus2){
  if(substr(elem, 1, 1)=="X"){cat(paste("\n",substr(elem,2,nchar(elem))))} 
  else{cat(paste("\n",elem))}}
cat("\n")



