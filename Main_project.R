######Packages######
#BiocManager packages installation for Volcano Plot#
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

#install development version
install.packages("devtools")
devtools::install_github("kevinblighe/EnhancedVolcano")

#Load the package into R session
library(EnhancedVolcano)
################
#string packages installation for using "str_replace_all" function for Sample Names of CNVs
install.packages("stringr")

library("stringr")
###############
# We use the "glmnet" package to perform feature selection for regression.
install.packages("glmnet")

library(glmnet) 
###############
##carData packages installation for using "avPlots" function for plotting Model Linear Regression
install.packages("carData")

###############################################################################
####Get Path####
DataPath = "E:/Code R"

####1.LUSC
path.LuscHealthy <- paste(DataPath,"lusc-rsem-fpkm-tcga_paired.txt",sep="/")
path.LuscCancer <- paste(DataPath,"lusc-rsem-fpkm-tcga-t_paired.txt", sep="/")
path.LuscCNV <- paste(DataPath,"lusc_CNV_core.txt",sep="/")


####2.Kidney
path.KidneyHealthy <- paste(DataPath,"kirc-rsem-fpkm-tcga_paired.txt",sep="/")
path.KidneyCancer <-paste(DataPath,"kirc-rsem-fpkm-tcga-t_paired.txt", sep="/")
path.KidneyCNV <- paste(DataPath,"kirc_CNV_core.txt",sep="/")

###############################################################################
source("Prepare_Data_Function.R")
source("Main_tests.R")
source("Regression_test.R")
######Data Preparation######
####1.LUSC
Data  <- Prepare_GE_Data(path.LuscHealthy,path.LuscCancer, path.LuscCNV)
Data_Lusc_Healthy = Data[[1]] #Gene Expression Healthy Data
Data_Lusc_Cancer = Data[[2]] #Gene Expression Cancer Data
Data_Lusc_CNV = Data[[3]]  #Copy Number Variations  Data
DC_Lusc_lm =Data[[4]] #intersection subset of Samples between CNVs and Cancer Data of GE 

#1.I. Hypothesis testing 
#1.I.a. paired testing , fold Change testing  and Volcano plotting
run_enhanced_volcano(Data_Lusc_Healthy, Data_Lusc_Cancer, "Lusc")

#1.I.b. independent Testing 
run_indep(Data_Lusc_Healthy, Data_Lusc_Cancer, "Lusc")

#Check Function
#DEGS and non-DEGS Comparison and output data summary 
run_intersect("Lusc", TRUE)

#1.I. Linear Regression
run_reg(DC_Lusc_lm,Data_Lusc_CNV,DataPath ,"Lusc" )
################################################################################
####2.Kidney
Data <- Prepare_GE_Data(path.KidneyHealthy,path.KidneyCancer, path.KidneyCNV)
Data_Kidney_Healthy = Data[[1]] #Gene Expression Healthy Data
Data_Kidney_Cancer = Data[[2]]  #Gene Expression Cancer Data
Data_Kidney_CNV = Data[[3]] #Copy Number Variations  Data
DC_Kidney_lm =Data[[4]]  #intersection subset of Samples between CNVs and Cancer Data of GE 

#2.I.Hypothesis testing 
#2.I.a. paired testing , fold Change testing  and Volcano plotting
run_enhanced_volcano(Data_Kidney_Healthy, Data_Kidney_Cancer, "Kidney")

#2.I.b. independent Testing 
run_indep(Data_Kidney_Healthy,Data_Kidney_Cancer, "Kidney")

#Check Function
#DEGS and non-DEGS Comparison and output data summary 
run_intersect("Kidney", TRUE)

#2.II. Linear Regression
run_reg(DC_Kidney_lm,Data_Kidney_CNV,DataPath , "Kidney")

### ### ### ### ### END of Program ### ### ### ### ###
