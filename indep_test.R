#● Do the samples follow the normal distribution?
#– Yes: Use the t-test
      # •So, we need to test for variance equality:
           # • If equal: use the t-test

indep_t_test_equal <- function(Data_Lusc_Cancer1,Data_Lusc_Healty1,GN_Lusc, cancer_type){
  
  p_value<- t.test(x= Data_Lusc_Cancer1, y=Data_Lusc_Healty1, paired = FALSE, alternative ='two.sided', var.equal = TRUE)$p.value
  
  if(p_value> 0.05 )#not DEGS
  {
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Indep_Not_DEGs_Report.txt",sep="_"),append=TRUE)
 
   } else #DEGS
  {
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Indep_DEGs_Report.txt",sep="_"),append=TRUE)
  }
}          # • If not equal: use the Wilch test.

indep_t_test_Not_equal <- function(Data_Lusc_Cancer1,Data_Lusc_Healty1,GN_Lusc, cancer_type){
 
   p_value<- t.test(x= Data_Lusc_Cancer1, y=Data_Lusc_Healty1, paired = FALSE, alternative ='two.sided', var.equal = FALSE)$p.value
  
  if(p_value> 0.05 )#not DEGS
  {
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Indep_Not_DEGs_Report.txt",sep="_"),append=TRUE)

  } else #DEGS
  {
       write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Indep_DEGs_Report.txt",sep="_"),append=TRUE)
  }
}


#– No: Use the Wilcoxon rank sum test, not Wilcoxon signed rank test.

indep_t_test_wilcox <- function(Data_Lusc_Cancer1,Data_Lusc_Healty1,GN_Lusc, cancer_type){
  
  p_value<- wilcox.test(x = Data_Lusc_Cancer1, y = Data_Lusc_Healty1, alternative = 'two.sided')$p.value
  
  if(p_value> 0.05 )#not DEGS
  {
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Indep_Not_DEGs_Report.txt",sep="_"),append=TRUE)

  } else #DEGS
  {
    write(paste(GN_Lusc,p_value, sep="\t"),file=paste(cancer_type,"Indep_DEGs_Report.txt",sep="_"),append=TRUE)
  }
}
