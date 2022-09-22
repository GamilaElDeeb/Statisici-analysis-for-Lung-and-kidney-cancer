#############################################################
#### Method_1. Hypothesis Testing ####
#For Case:1. Samples are paired.
#Data Distribution Testing
run_paired<- function(Data_Healthy, Data_Cancer, cancer_type){
    print(paste("Run Paired for: ",cancer_type))

    ##To check for normality for the paired difference
    paired.Difference <- Data_Healthy - Data_Cancer 
    Shapiro.Pvalue <- apply(paired.Difference,1,function(x) shapiro.test(x)$p.value)
    
    source("paired_test.R")
    
    allnormal <- 1
    for(i in 1:length(Shapiro.Pvalue)){
        if(Shapiro.Pvalue[i] < 0.05){
            allnormal <- 0
            break
        }
    }
    if(allnormal == 1){
        print("Data Normally Distributed, Applying T_Test ...")
        p_value <- apply_ttest(Data_Cancer, Data_Healthy, Shapiro.Pvalue, cancer_type)
    }else{
        print("Data Not-Normally Distributed, Applying Wilcox_Test ...")
        p_value <- apply_wilcox(Data_Cancer, Data_Healthy, Shapiro.Pvalue, cancer_type)
    }
    print(paste("Paired Ended for:",cancer_type))
    
    return(p_value)
}


##For Case:2. Samples are independent.
#Data Distribution Testing

run_indep<- function(Data_Healthy, Data_Cancer, cancer_type){
    print(paste("Run Indep for: ",cancer_type))

    source("indep_test.R")
    
    ##To check for normality for the GE Data
    Shapiro.indepH.Pvalue <- apply(Data_Healthy,1,function(x) shapiro.test(x)$p.value)
    Shapiro.indepD.Pvalue <- apply(Data_Cancer,1,function(x) shapiro.test(x)$p.value)

    indep.allNormal <- 1
    for(i in 1:length(Shapiro.indepD.Pvalue)){
        if(Shapiro.indepH.Pvalue[i] < 0.05 || Shapiro.indepD.Pvalue[i] < 0.05){
            indep.allNormal <- 0
            break
        }
    }
    if(indep.allNormal == 1){
        print("Indep: Data Normally Distributed..")
        equalPValue <- 1
        for(i in 1:nrow(Data_Healthy)){
            Data_Cancer1 <- as.numeric(Data_Cancer[i,c(1:50)])
            Data_Healthy1 <- as.numeric(Data_Healthy[i,c(1:50)])
            p1_value <- var.test(Data_Cancer1,Data_Healthy1)$pvalue
            if(p1_value > 0.05){
                equalPValue <- 0
                break
            }
        }
        #for normal distribution , variance check..
        if(equalPValue == 1){
            print("Indep: All Normal - Equal Variance , Applying T_Test ...")
            for(i in 1:nrow(Data_Healthy)){
                Data_Cancer1 <- as.numeric(Data_Cancer[i,c(1:50)])
                Data_Healthy1 <- as.numeric(Data_Healthy[i,c(1:50)])
                GN_Names <- row.names(Data_Cancer[i,])
                indep_t_test_equal(Data_Cancer1,Data_Healthy1,GN_Names, cancer_type)
            }
        }else{
            print("Indep: All Normal - Variance Not Equal , Applying Wilch test ...")
            for(i in 1:nrow(Data_Healthy)){
                Data_Cancer1 <- as.numeric(Data_Cancer[i,c(1:50)])
                Data_Healthy1 <- as.numeric(Data_Healthy[i,c(1:50)])
                GN_Names <- row.names(Data_Cancer[i,])
                indep_t_test_Not_equal(Data_Cancer1,Data_Healthy1,GN_Names, cancer_type)
            }
        }
    }else{
        print("Indep: Data Not-Normally Distributed , Appling Wilcoxon rank sum test ...")
        for(i in 1:nrow(Data_Healthy)){
                Data_Cancer1 <- as.numeric(Data_Cancer[i,c(1:50)])
                Data_Healthy1 <- as.numeric(Data_Healthy[i,c(1:50)])
                GN_Names <- row.names(Data_Cancer[i,])
                indep_t_test_wilcox(Data_Cancer1,Data_Healthy1,GN_Names, cancer_type)
            }
    }
    print(paste("Indep Ended for:",cancer_type))
}   

##DEGS and Not-DEGS for Independent and paired 
##Report how different these two sets of genes :

run_intersect <- function(cancer_type, print_check){
  print(paste("Intersection for: ",cancer_type))
   
  #Data Path Determination and Reading Files For DEGS
  Indep.path <- paste(cancer_type,"Indep_DEGs_Report.txt",sep="_")
  Paired.path <- paste(cancer_type,"Paired_DEGs_Report.txt",sep="_")
  
  Indep.df <- read.delim(Indep.path, row.names = 1, header=FALSE)
  Paired.df <- read.delim(Paired.path, row.names = 1, header=FALSE)
  
  Indep_rowname <- rownames(Indep.df)
  Paired_rowname <- row.names(Paired.df)
  
  #Extraction Common  the differentially expressed genes
  Common <- intersect(Indep_rowname,Paired_rowname)
  lapply(Common, write, paste(cancer_type,"Common_DEGs.txt",sep="_"), append=TRUE)
  
  #Extraction of Unique DEGS For Independent Data
  diffIndpOnly <- setdiff(Indep_rowname,Common)
  lapply(diffIndpOnly, write, paste(cancer_type,"Unique_DEGS_Independent.txt",sep="_"), append=TRUE)
  
  #Extraction of Unique DEGS For Paired Data
  diffPairedOnly <- setdiff(Paired_rowname,Common)
  lapply(diffPairedOnly, write, paste(cancer_type,"Unique_DEGS_Paired.txt",sep="_"), append=TRUE)
  
  #Check Point and Data Summary
  if(print_check){
    print(paste("Indep has Genes : "),length(Indep_rowname))
    
    print(paste("Paired has Genes : ",length(Paired_rowname)))

    print(paste("For Common Genes Length :",length(Common)))

    print(paste("For Unique Genes Indep Length : ",length(diffIndpOnly)))

    print(paste("For Unique Genes Paired Length : ",length(diffPairedOnly)))
  }
  
}

#############################################################
#### Method_2. Fold change ####

run_fc <- function(Data_Healthy, Data_Cancer, cancer_type){
    print(paste("Starting Fold change for : ",cancer_type))

    GN <- as.character(rownames(Data_Healthy))
    
    #Fold change Data Preparation
    healthyMean <- rowMeans(Data_Healthy)
    cancerMean <- rowMeans(Data_Cancer)
    fc <- cancerMean/healthyMean
    
    fclog2 <- log2(fc)
    fclog2abs <- abs(fclog2)
    threshold <- log2(1.5)
    
    #Extraction DEGS and Not-DEGS Using Fold change
    for(i in 1:length(fclog2abs)){
        if(fclog2abs[i]>= threshold){
            write(paste(fclog2[i],GN[i], sep="\t"),file=paste(cancer_type,"FC_DEGs_Report.txt",sep="_"),append=TRUE)
        }else{
            write(paste(fclog2[i],GN[i], sep="\t"),file=paste(cancer_type,"FC_Not_DEGs_Report.txt",sep="_"),append=TRUE)
        }
    }
    print("Ending  Fold change .")
    
    return(fclog2)
    
}

#############################################################
#### Method_3. Both of them (volcano plot) ####

run_enhanced_volcano <- function(Data_Healthy, Data_Cancer ,cancer_type){
    
    print(paste("Starting Volcano Plot for",cancer_type))
  
    #volcano plot Data Preparation
    pvalue <- run_paired(Data_Healthy, Data_Cancer, cancer_type)
    log2FoldChange <- run_fc(Data_Healthy, Data_Cancer, cancer_type)
    res <- data.frame( log2FoldChange,pvalue)
    
    # The default cut-off for log2FC is >|2|; the default cut-off for P value is 10e-6.

    EnhancedVolcano(res,
                    lab = rownames(res),
                    x = 'log2FoldChange',
                    y = 'pvalue',
                    pCutoff = 0.05,   #alpha - on Y_axis
                    FCcutoff = 0.6    #Threshold - on X_axis
    )
    
    #Extraction  Names of DEGS Using volcano plot
    out <- EnhancedVolcano(res,
                           lab = rownames(res),
                           x = 'log2FoldChange',
                           y = 'pvalue',
                           pCutoff = 0.05,   #alpha - on Y_axis
                           FCcutoff = 0.6,    #Threshold - on X_axis
    )
    Volcano_DEGS <- out$layers[[4]]$data$lab
    
    #Preparation DEGS For GSEA Software
    Volcano_DEGS_DF_Healthy <- Data_Healthy[rownames(Data_Healthy) %in% Volcano_DEGS, ] 
    Volcano_DEGS_DF_Cancer <- Data_Cancer[rownames(Data_Cancer) %in% Volcano_DEGS, ]
    DEGS_ChangeN_Cancer <- sub("T", "H", colnames(Volcano_DEGS_DF_Cancer)) 
    colnames(Volcano_DEGS_DF_Cancer) <-DEGS_ChangeN_Cancer
    Volcano_DEGS_DH <-data.frame( Volcano_DEGS_DF_Cancer,Volcano_DEGS_DF_Healthy )
    write.table(Volcano_DEGS_DH, paste(cancer_type,"Volcano_DEGs.txt",sep="_"),row.names = TRUE,sep = "\t")
    

    print("Ending  Volcano Plot .")
}
### ### ### ### ### END of Requirement.1 ### ### ### ### ###
