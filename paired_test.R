# ● Does the paired sample follow the normal distribution?
# – Yes: Use the t-test

apply_ttest <- function(Data_Cancer, Data_Healty, Shapiro.Pvalue, cancer_type){
  
  all_pvalues <- c()
  for(i in 1:nrow(Data_Cancer)){
    
    #convert from data frame as numeric vector
    Data_Cancer1 <- as.numeric(Data_Cancer[i,c(1:50)])
    Data_Healty1 <- as.numeric(Data_Healty[i,c(1:50)])
    
    GN_Lusc <- row.names(Data_Cancer[i,])
    
    Shapiro.pvalue <- Shapiro.Pvalue.Lung[i]
    
    p1_value <- t_test(Data_Cancer1,Data_Healty1,GN_Lusc, cancer_type)
    
    all_pvalues = c(all_pvalues, p1_value)        
  }
  return(all_pvalues)
}

t_test =  function(Data_Cancer1,Data_Healty1,GN_Lusc, cancer_type)
{
  p_value<- t.test(Data_Cancer1,Data_Healty1, paired = TRUE, alternative = "two.sided")$p.value
  
  if(p_value> 0.05 )#not DEGS
  {
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Paired_Not_DEGs_Report.txt",sep="_"),append=TRUE)
  } else #DEGS
  {
    p_value<- t.test(Data_Cancer1,Data_Healty1, paired = TRUE, alternative = "two.sided")$
    
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Paired_DEGs_Report.txt",sep="_"),append=TRUE)
  }
  return(p_value)
}


# – No: Use the Wilcoxon signed rank test.

apply_wilcox <- function(Data_Cancer,Data_Healty, Shapiro.Pvalue, cancer_type){
 
   all_pvalues <- c()
  for(i in 1:nrow(Data_Cancer)){
  
    Data_Cancer1 <- as.numeric(Data_Cancer[i,c(1:50)])
    Data_Healty1 <- as.numeric(Data_Healty[i,c(1:50)])
    
    GN_Lusc <- row.names(Data_Cancer[i,])
    
    Shapiro.pvalue <- Shapiro.Pvalue[i]
    
    p1_value <- wilcox_test(Data_Cancer1,Data_Healty1,GN_Lusc, cancer_type)
    
    all_pvalues = c(all_pvalues, p1_value)    
  }
   
  return(all_pvalues)
}

wilcox_test =  function(Data_Cancer1,Data_Healty1,GN_Lusc, cancer_type)
{
  p_value<- wilcox.test(Data_Cancer1,Data_Healty1, paired = TRUE, alternative = "two.sided")$p.value
  
  if(p_value> 0.05 )#not DEGS
  {
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Paired_Not_DEGs_Report.txt",sep="_"),append=TRUE)
 
   } else #DEGS
  {
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Paired_DEGs_Report.txt",sep="_"),append=TRUE)
  }
  
  return(p_value)
}
### ### ### ### ### END of Paired ### ### ### ### ###

