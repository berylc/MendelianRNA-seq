library(plyr)
library(dplyr)
source("../QC/bin.R")

#Read in GTEx expression matrix, attributes and phenotype files. 
gtex_rpkm <- read.delim("../data/gtex_source_data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct.gz",stringsAsFactors=F, skip=2)
gtex_attributes <- read.delim("../data/gtex_source_data/GTEx_Data_V6_Annotations_SampleAttributesDS.txt", stringsAsFactors = F) 
gtex_attributes$SAMPID <- gsub("-","\\.",gtex_attributes$SAMPID )
gtex_phenotypes <- read.delim("../data/gtex_source_data/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt",stringsAsFactors=F)[,1:2]


#Note that resulting files from this snippet are already available under MendelianRNA-seq/data/ so you don't need to run this code to work through the example. This snippet is for those wishing to make their own tissue check files. Just make sure to change any relevant file paths in the QC code with the resulting files here.  

#This section is to make files required for MuscleCheck.R. 
#Extract tissues of interest
tissue_attributes<- gtex_attributes %>% filter(SMTSD %in% c("Adipose - Subcutaneous", "Muscle - Skeletal","Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)","Cells - Transformed fibroblasts")) %>% select(SAMPID,SMTSD)

tissue_attributes$SMTSD <- gsub(" ","",tissue_attributes$SMTSD )

tissue_samples <- intersect(tissue_attributes$SAMPID,names(gtex_rpkm))
tissue_rpkm<-gtex_rpkm[,c("Name","Description",tissue_samples)]

#Subset expression matrix to only genes that are tissue preferentially expressed. If you change the tissue selection, make sure to change the tissue_pref_genes to the relevant genes from the Mele et al supplementary data. 
tissue_rpkm<- extractGenesOfInterest(tissue_rpkm, tissue_pref_genes)

#Rename samples to reflect tissues
tissue_rpkm <-merge(tissue_rpkm, tissue_attributes, by.x='row.names',by.y='SAMPID', all.x=T, all.y=F)
tissue_rpkm$SMTSD<- gsub(")$","",tissue_rpkm$SMTSD)
rownames(tissue_rpkm) <- make.names(tissue_rpkm$SMTSD, unique=T)
tissue_rpkm<-subset(tissue_rpkm, select = -c(Row.names,SMTSD))

write.table(tissue_rpkm,"./my.own.gtex_expression_tissue_preferential_fibs.msck.skn.adbp.txt",quote=F, col.names=T, row.names=T, sep="\t")

#This section is to make files required for SexCheck.R 

#Since we are comparing to muscle, muscle samples from the expression matrix 
muscle_samples <- subset(gtex_attributes, SMTSD == "Muscle - Skeletal")$SAMPID
muscle_samples <- intersect(muscle_samples, names(gtex_rpkm))
sex_rpkm<-gtex_rpkm[,c("Name","Description",muscle_samples)]

#Subset expression matrix to only genes that are sex preferentially expressed.
sex_rpkm<- extractGenesOfInterest(sex_rpkm, sex_biased_genes$gene_id)

#Annotate the GTEx file with relevant sex information of samples
rownames(sex_rpkm) <- sapply(strsplit(rownames(sex_rpkm), split="\\."), function(x) paste(x[1],x[2],sep="-"))
sex_rpkm<-merge(sex_rpkm, gtex_phenotypes, by.x='row.names', by.y = 'SUBJID', all.x=T, all.y=F)
sex_rpkm$GENDER <- recode(sex_rpkm$GENDER, `1`="male", `2`="female")
sex_rpkm$Row.names <- make.names(sex_rpkm$GENDER, unique=T)
sex_rpkm<- subset(sex_rpkm, select = -GENDER)

write.table(tissue_rpkm,"./my.own.gtex_expression_sex_biased_genes.txt",quote=F, col.names=T, row.names=T, sep="\t")
