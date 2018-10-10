## ---- eval=FALSE---------------------------------------------------------
#  
#  library(miRLAB)
#  
#  dataset=system.file("extdata", "EMT35.csv", package="miRLAB")
#  cause=1:35 #column 1:35 are miRNAs
#  effect=36:1189 #column 36:1189 are mRNAs
#  
#  #predict miRNA targets using Pearson correlation
#  pearson=Pearson(dataset, cause, effect)
#  
#  #predict miRNA targets using Mutual Information
#  mi=MI(dataset, cause, effect)
#  
#  #predict miRNA targets using causal inference
#  ida=IDA(dataset, cause, effect, "stable", 0.01)
#  
#  #predict miRNA targets using linear regression
#  lasso=Lasso(dataset, cause, effect)
#  

## ---- eval=TRUE----------------------------------------------------------
library(miRLAB)
#validate the results of the top100 targets of each miRNA predicted 
#by the four methods
dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
pearson=Pearson(dataset, 1:3, 4:18)
miR200aTop10=bRank(pearson, 3, 10, TRUE)
groundtruth=system.file("extdata", "Toygroundtruth.csv", package="miRLAB")
miR200aTop10Confirmed = Validation(miR200aTop10, groundtruth)

## ---- eval=TRUE----------------------------------------------------------
library(miRLAB)
#validate the results of the top100 targets of each miRNA predicted 
#by the four methods
dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
EMTresults=Pearson(dataset, 1:3, 4:18)
top10=Extopk(EMTresults, 10)
groundtruth=system.file("extdata", "Toygroundtruth.csv", package="miRLAB")
top10Confirmed = Validation(top10, groundtruth)

## ---- eval=TRUE----------------------------------------------------------
library(miRLAB)
dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
dataset=Read(dataset)
dataset[1:5,1:7]

## ---- eval=TRUE----------------------------------------------------------
library(miRLAB)
groundtruth=system.file("extdata", "Toygroundtruth.csv", package="miRLAB")
groundtruth=Read(groundtruth)
groundtruth[1:5,]

