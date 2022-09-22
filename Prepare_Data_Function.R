Prepare_GE_Data <- function(path.Healthy , path.Cancer, path.CNV) {
  
  ####Read GE Data in any extension using read.delim & set Row names column######
  GE.Healthy<- read.delim(path.Healthy, row.names = 1, header=TRUE)
  GE.Cancer <- read.delim(path.Cancer, row.names = 1, header=TRUE)
  
  ###Read CNV Data
  CNV.Core <- read.delim(path.CNV, row.names = 1, header=TRUE)

  #Data Filtration for GE and CNVs
  #####Exclude Genes of more than 50% zeros using Median Average######
  Filter.DataHealthy<-GE.Healthy[which(apply(GE.Healthy, 1, median)!=0),]
  Filter.DataCancer<-GE.Cancer[which(apply(GE.Cancer, 1, median)!=0),]
  Filter.CNV<-CNV.Core[ ,which(apply(CNV.Core, 2, median)!=0)]
  
  ###Remove the Gene ID Column in GE#####
  Filter.DataHealthy<-Filter.DataHealthy[,-grep("Entrez_Gene_Id",colnames(Filter.DataHealthy))]
  Filter.DataCancer<-Filter.DataCancer[,-grep("Entrez_Gene_Id",colnames(Filter.DataCancer))]
  
  ###Getting Common Row names After Filtration in GE Data#######
  Healthy_rowname <- rownames(Filter.DataHealthy)
  Cancer_rowname <- row.names(Filter.DataCancer)
  Common <- intersect(Healthy_rowname,Cancer_rowname)
  ######Extraction of common GE data rows after Filtration
  Data_subset_Healty <- Filter.DataHealthy[rownames(Filter.DataHealthy) %in% Common, ]  
  Data_subset_Cancer <- Filter.DataCancer[rownames(Filter.DataCancer) %in% Common, ] 
  
  #############
  #Data Preparation and Filtration  of samples for linear regression
  #row name modification in CNV
  CNV.F_rowname <- row.names(Filter.CNV) 
  CNV.F.replace_rowname<- str_replace_all(CNV.F_rowname,"-",".")
  row.names(Filter.CNV) <- CNV.F.replace_rowname
  CNV.F_Arowname <- row.names(Filter.CNV)
  #Getting Common Samples between cancer GE and CNVs
  Cancer_colname <-colnames(Data_subset_Cancer)
  Common.GE_C.CNV <- intersect(CNV.F_Arowname,Cancer_colname)
  Data_CNV <- Filter.CNV[rownames(Filter.CNV) %in% Common.GE_C.CNV, ]
  Data_Cancer_lm <-Data_subset_Cancer[ ,colnames(Data_subset_Cancer) %in% Common.GE_C.CNV] 
  
  return(list(Data_subset_Healty,Data_subset_Cancer,Data_CNV,Data_Cancer_lm))
  
}