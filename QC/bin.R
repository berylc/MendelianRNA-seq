tissue_pref_genes<-read.delim("/humgen/atgu1/fs03/berylc/MuscDisease/Git/MendelianRNA-seq/data/TissuePreferentiallyExprGenes_FB_MSCK_SKN_ADBP.txt",header=T,stringsAsFactors=F)$gene

extractGenesOfInterest<-function(rpkm_file,makecolnames=T,gene_list,patient_col=3,type=NA){
  rpkm_file<-subset(rpkm_file, Description %in% gene_list)
  rpkm_file$Name<-NULL
  rownames(rpkm_file)<-rpkm_file$Description
  rpkm_file$Description<-NULL
  
  if(makecolnames==T){colnames(rpkm_file)<-make.names(rep(type,ncol(rpkm_file)),unique=T)}
  rpkm_file <- data.frame(t(rpkm_file),stringsAsFactors=F)
  return(rpkm_file)}



getSexBiased_xist_y<- function(rpkmfile,patientcolbegin=3,out,sexBiasedGenes_yChrom){
  XIST<-rpkmfile[rpkmfile$Description=="XIST",patientcolbegin:ncol(rpkmfile)]
  gtexYchrom<-data.frame(rpkmfile[rpkmfile$Description %in% sexBiasedGenes_yChrom$gene_name,patientcolbegin:ncol(rpkmfile)])
  add_row<-nrow(gtexYchrom)+1
  gtexYchrom[add_row,]<-apply(gtexYchrom,2,mean)
  gtexYchromVals<-gtexYchrom[add_row,]
  XIST<-data.frame(t(XIST),stringsAsFactors = F)
  gtexYchromVals<-data.frame(t(gtexYchromVals),stringsAsFactors = F)
  both<<-merge(XIST,gtexYchromVals,by="row.names")
  names(both)<-c("sample","XIST","avgYchrom")
  assign(x=paste("XIST_Y",out,sep="_"),value=both,envir=globalenv())
}

getSexBiasedAll <- function(rpkmfile,patientcolbegin=3,out,sexBiasedGenes){
  SexBiased_all<-rpkmfile[rpkmfile$Description %in% sexBiasedGenes$gene_name,]
  assign(x=paste("AllSexBiasedGenes",out,sep=""),value=SexBiased_all,envir=globalenv())
}

getcombinedPCAData<- function(rpkmfilegtex,rpkmfilepatients,out,sexBiasedGenes){   
  SexBiased_gtex<-rpkmfilegtex[rpkmfilegtex$Description %in% sexBiasedGenes$gene_name,2:ncol(rpkmfilegtex)]     
  SexBiased_patients<-rpkmfilepatients[rpkmfilepatients$Description %in% sexBiasedGenes$gene_name,2:ncol(rpkmfilepatients)]   
  SexBiased_all<-merge(SexBiased_patients,SexBiased_gtex,by="Description")   
  SexBiased_all<-SexBiased_all[!duplicated(SexBiased_all$Description),]   
  SexBiased_all<-SexBiased_all[,2:ncol(SexBiased_all)]   
  SexBiased_all <- SexBiased_all + 1   
  SexBiased_all <- log2(SexBiased_all)   
  SexBiased_all<-data.frame(t(SexBiased_all),stringsAsFactors = F)  
  thrTissPCA<-prcomp(na.omit(SexBiased_all), retx=TRUE,na.action=na.omit,center=TRUE)   
  assign(x=paste("AllSexBiased_RPKM",out,sep="_"),value=SexBiased_all,envir=globalenv()) 
  assign(x=paste("PCADat",out,sep="_"),value=thrTissPCA,envir=globalenv()) }

#After initial precompution, I changed this to work with patient files. Need to add Klinefelters samples to have it work again. 
addPhenotype<- function (All,out,gtex_phenotypes){
  PhenIds<-gtex_phenotypes$SUBJID
  firstCol<-ncol(All)+1
  secCol<-ncol(All)+2
  for(m in 1:length(rownames(All))){
    mS = strsplit(as.character(rownames(All)[m]),split="\\.")[[1]]
    newID = paste(mS[1],"-",mS[2],sep="")
    All[m,firstCol]<-newID
    if(newID %in% PhenIds){
      GenderCode=gtex_phenotypes[gtex_phenotypes$SUBJID==newID,"GENDER"]
      if(GenderCode==1){All[m,secCol]<-"Male"}
      if(GenderCode==2){All[m,secCol]<-"Female"}
    }
  #  if(rownames(All)[m] %in% KlinefeltersSamples){
   #   All[m,secCol]<-"Klinefelter's"
  #  }
    if(!newID %in% PhenIds){All[m,secCol]<-"Patient"}
    
  }
#  names(All)<-c("sample","XIST","avgYchrom","shortID","gender")
  assign(x=paste("WithSex",out,sep="_"),value=All,envir=globalenv())
}


