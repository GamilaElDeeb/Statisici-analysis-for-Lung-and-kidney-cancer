run_Feature_Selection <- function(Data.CNV , One.DEGS){
  fit.cv <- cv.glmnet(Data.CNV, One.DEGS, family="gaussian", alpha=1, standardize=FALSE, nfolds=5)
  # Getting the value of lowest value of lambda.
  lambda <- fit.cv$lambda.min  
  # Here we compute the variables regression coefficeints after being penalized (set to zero).
  model <- glmnet(Data.CNV, One.DEGS, family="gaussian", alpha=1, lambda=lambda, standardize=FALSE)
  # exclude intercept (the first value).
  coef.fit <- coef(model, s=lambda)[2:(ncol(Data.CNV)+1)]   
  # Keep features with non-zero regression coefficients.
  features.in <- which(abs(coef.fit) > 0) 
  # The variables matrix after removng the penalized featurs.
  FS.CNV = Data.CNV[,features.in]  
  
  return(FS.CNV)
}

run_reg <- function(Data_Cancer_lm,Data_CNV,DataPath,cancer_type){
    
    path.DEGS.Paired <- paste(DataPath,paste(cancer_type,"Paired_DEGs_Report.txt",sep="_"),sep="/")

    read.DEGS.Paired<- read.delim(path.DEGS.Paired,row.names=NULL, header=FALSE)

    colnames(read.DEGS.Paired) <- c("DEGS","Pvalue")
     #to obtain the name of DEGS
    Min.5_DEGS <- head(read.DEGS.Paired[order(read.DEGS.Paired$Pvalue),],5)$DEGS
    print(paste("Name of Five DEGS: ", Min.5_DEGS )) 
     #to get Samples of this DEGS
    Min.5_DEGS_Sample <- as.matrix(Data_Cancer_lm[Min.5_DEGS, ])
    Data_CNV <- as.matrix(Data_CNV)
    
    for(i in 1:nrow(Min.5_DEGS_Sample)){
      
      if(ncol(Data_CNV) > nrow(Data_CNV)){
        
        #Feature Selection --Lung
        Data_CNV <-run_Feature_Selection(Data_CNV ,Min.5_DEGS_Sample[i,]) 
      }
        DEGS.1<- Min.5_DEGS_Sample[i,]
        model <- lm(DEGS.1 ~ Data_CNV) # Linear regression.
        print(summary(model))
        capture.output(summary(model), file = "Kidney_Regression_Report.txt", append=T)
        require(car)
        avPlots(model)

      
    }
}


