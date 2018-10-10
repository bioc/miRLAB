# This file contains a collection of methods for inferring miRNA-mRNA regulatory relationships
# Each method will take the dataset (.csv) as an input with the first n colomns are miRNAs 
# and the last m columns are mRNAs. Rows are samples. We need to specify the number of miRNAs in the dataset
# and the index of the miRNA, i.e. column that contains the miNRA. Topk is the number of top k interactions we
# would like to extract.






#getData from TCGA by TCGAbiolinks package
#' @import TCGAbiolinks
#' @import dplyr
#' @import SummarizedExperiment
#' @importFrom utils read.csv
#' @importFrom utils write.csv
#' @importFrom methods new
#' @importFrom stats p.adjust
#' @importFrom stats cor
#' @importFrom stats cov
#' @importFrom stats rnorm
#' @importFrom stats cancor
#' @importFrom stats phyper
#' @importFrom stats median
#' @import RCurl
#' @import httr
#' @import stringr
#' @import ctc
#' @import heatmap.plus


`%:::%` = function(pkg, fun) get(fun, envir = asNamespace(pkg), inherits = FALSE)

# A list of normal/tumor sample types, you can get 19 types via the command: TCGAbiolinks:::getBarcodeDefinition()$tissue.definition

TCGAdownload = function(CancerType="BRCA", Datatype="miRNA", Platform="miRNA-Seq", SampleType="tumor", ncases) {
  normal_pos = c(10:14)
  # tumor_barcode = TCGAbiolinks:::getBarcodeDefinition()$tissue.definition[-normal_pos]
  # normal_barcode = TCGAbiolinks:::getBarcodeDefinition()$tissue.definition[normal_pos]
  x <- 'TCGAbiolinks' %:::% 'getBarcodeDefinition'
  tumor_barcode = x()$tissue.definition[-normal_pos]
  normal_barcode = x()$tissue.definition[normal_pos]
  
  #find the projects conducted by The Cancer Genome Atlas genome program 
  project = paste("TCGA-", CancerType, sep = "")
  
  sample_type = switch(tolower(SampleType), 
                       "tumor" = tumor_barcode,
                       "normal" = normal_barcode)
  
  switch(tolower(Datatype), 
         "mirna" = { 
           #data_category = "Transcriptome Profiling"
           data_category = "Gene expression"
           data_type = "miRNA gene quantification"
           #data_type="miRNA Expression Quantification"
           file_type = "mirbase20.mirna" 
         },
         
         "mrna" = { 
           #data_category = "Gene expression"
           data_category="Gene expression"
           data_type = "Gene expression quantification"
           #file_type = ".rsem.genes.normalized_results"
           file_type = "normalized_results"
         },
         
         "somatic mutation" = { 
           data_category = "Simple nucleotide variation"
           data_type = "Simple somatic mutation"
           file_type = ".automated.somatic.maf"
         },
         
         "clinical" = {
           result = GDCquery_clinic(project)
           return(result)
         },
         
         "methylation" = { #get data of methylation
           data_category = "DNA methylation"
           data_type = "Methylation beta value"
           file_type = FALSE
         }
  )
  
  #eg. function argument Platform = RNA-Seq = experimental.strategy of GDCdownload
  experimental_strategy = Platform
  
  #query all of the cases based on the project and their types of data
  query = GDCquery(project = project, 
                   data.category =data_category,
                   data.type =data_type,
                   platform = "Illumina HiSeq",
                   #platform = "IlluminaGA_RNASeqV2",
                   file.type = file_type,
                   legacy = TRUE, access = "open", 
                   experimental.strategy = experimental_strategy,
                   #sample.type = sample_type,
                   #barcode=c("TCGA-A2-A1FZ-01A-51R-A14D-07")
  )
  print(query$result[[1]])
  # get barcode for the first [ncases] cases and download them, if ncases = NA then download all
  if (!missing(ncases)) {
    getcases = getResults(query)$cases[1:ncases]
    
    query = GDCquery(project = project,
                     data.category = data_category,
                     data.type = data_type,
                     file.type = file_type,
                     legacy = TRUE, access = "open",
                     experimental.strategy = experimental_strategy,
                     sample.type = sample_type,
                     barcode = getcases
    )
  }
  
  GDCdownload(query)
  
  result_tmp =GDCprepare(query)
  #output
  switch(tolower(Datatype), 
         "mirna" = { 
           rownames(result_tmp) = t(as.character(result_tmp[,1]))
           result_tmp = as.matrix(result_tmp[,grep("reads_per_million", colnames(result_tmp))])
           colnames(result_tmp) = substring(colnames(result_tmp), 32)
           return(result_tmp)
         },
         
         "mrna" = { 
           return(assay(result_tmp, "normalized_count")) 
         }, 
         
         "somatic mutation" = {
           return(result_tmp)
         },
         
         "methylation" = {
           return(assay(result_tmp))
         }
  )
}


#' getData from GDC
#' @param cancerName The name of cancer in string format
#' @return dataset in matrix format
#' @examples
#' RawGBM=getData("GBM")
#' @export
getData=function(cancerName){
  mRNA = TCGAdownload(CancerType=cancerName,Datatype = "mRNA", Platform = "RNA-Seq")
  miRNA = TCGAdownload(CancerType=cancerName,Datatype = "miRNA", Platform = "miRNA-Seq")
  unlink("GDCdata",recursive = TRUE)
  unlink("MANIFEST.txt",recursive = TRUE)

  samples_mRNA = colnames(mRNA)
  samples_miRNA = colnames(miRNA)
  match = match(substring(samples_miRNA, 1, 20), substring(samples_mRNA, 1, 20))
  
  getmatch_data = samples_mRNA[match]
  match_data_names = cbind(samples_miRNA, getmatch_data)
  
  ans = c() #ans is our goal of combination mRNA and miRNA
  
  for (i in c(1: length(getmatch_data))) {
    if (!is.na(match_data_names[i, 2])) {
      #sample_name = substring(match_data_names[i, 1], 1, 20)
      tmp = c(as.numeric(miRNA[, i]), as.numeric(mRNA[, match[i]])) #data of sample
      ans = rbind(ans, tmp)
    }
  }
  
  empty = which(is.na(match_data_names[, 2]))
  match_data_names = match_data_names[-empty, ]
  rownames(ans) = NULL
  colnames(ans) = c(rownames(miRNA), rownames(mRNA))
  return(ans)
}





#' Read dataset from csv file
#' @param dataset The input dataset in csv format
#' @return dataset in matrix format
#' @examples
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' data=Read(dataset)
#' @export
Read<-function(dataset){
		data<-read.csv(dataset, header=TRUE, sep=",")
		#data<-data[,-1]
		return(data)
	}

############ Get data header from dataset ###############
#' Read the header of the dataset
#' @return the header of the dataset
#' @param dataset the character string of the names of the dataset in csv format, e.g. "ToyEMT.csv"
#' @examples
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' header=readHeader(dataset)
#' @export
readHeader<-function(dataset){
  data<-read.csv(dataset, header=FALSE)
  header<-character()
  for (i in 1:ncol(data)){
  header[i]=toString(data[1,i])
  }
  return(header)
}

############ Convert target information into a binary matrix, the rows are mRNAs and the columns are miRNAs ###############
queryTargetFile<-function(miRNA,mRNA,file){
  #the column name of file should be tagged as "mir" and "gene"
  data<-Read(file)
  mir=as.character(data$mir)
  #   mir<-paste("hsa-",sep="",mir);mir<-sub('r','R',mir)
  gene=as.character(data$gene)
  #   symbol<-geneSymbol(gene)
  sum=0
  rep<-replicate(length(miRNA),mir)
  edge=matrix(FALSE,length(miRNA)+length(mRNA),length(miRNA)+length(mRNA))
  for (i in 1:length(mir)){
    #     print(i)
    if (length(which(rep[i,]==miRNA)>0)){
      match1<-which(rep[i,]==miRNA,arr.ind=TRUE)
      #gene: gene[i] #mirna: miRNA[match]
      rep2<-replicate(length(mRNA),gene[i])
      match2<-which(rep2==mRNA,arr.ind=TRUE)
      edge[match1,match2+length(miRNA)]=TRUE
    }
  }
  return(edge)
}
	
	
#' Stardarsise the dataset
#' Stadardise the dataset to have mean=0 and std=1 in each column.
#' @param dataset The input dataset in csv format. e.g. "ToyEMT.csv". The first column is the sample name.
#' @return The standardised dataset.
#' @examples
#' \dontrun{
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' stdata=Standardise(dataset)
#' }
	Standardise<-function(dataset){
		ncol<-ncol(dataset)
		nrow<-nrow(dataset)
		stdData<-matrix(nrow=nrow,ncol=ncol-1)
              rownames(stdData)<-dataset[,1]
              colnames<-colnames(dataset)
              colnames(stdData)<-colnames[2:ncol]
		for (i in 2:ncol){
			stdData[,i-1]<-scale(dataset[i],center=TRUE, scale=TRUE)
			}
		return(stdData)
	}

#' Filter, impute, and normalise data.
#' 
#' Remove the genes (rows) that have more than r\% of missing data;
#' use the impute package to fill in missing data, and finally normalise the data.
#' @importFrom impute impute.knn
#' @importFrom limma normalizeBetweenArrays
#' @param dataset The input dataset in csv format. e.g. "EMT.csv" 
#' @param r The rate threshold to filter the records (genes). Genes with more than r\% missing data will be removed.
#' @return The processed dataset.
#' @references 
#' 1. Hastie T, Tibshirani R, Narasimhan B and Chu G. impute: Imputation for microarray data. R package version 1.42.0.
#' 
#' 2. Smyth, G.K. (2005). Limma: linear models for microarray data. In Bioinformatics and computational biology solutions using R and Bioconductor (pp. 397-420). Springer New York.
#' @examples
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' impdata=ImputeNormData(dataset, 0.1)
#' @export
ImputeNormData<-function(dataset, r){
               #library(impute)
               Rawdata<-Read(dataset)
               #Rawnames<-as.matrix(Rawdata)[,1]
              # Rawdata<-Rawdata[,-1]
               #rownames(Rawdata)<-Rawnames
                              
               #Removing genes which samples having more r% zero counts or missing value
               Rawdata<-Rawdata[which(rowSums(Rawdata==0)<r*dim(Rawdata)[2]),]
               data<-Rawdata[which(rowSums(is.na(Rawdata))<r*dim(Rawdata)[2]),]
             
               #if(exists(".Random.seed")) rm(.Random.seed)
               imputed<-impute.knn(as.matrix(data) ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
               dataimputed<-imputed$data
               
               #Make log2 transformation for expression data
               dataimputed<-log2(dataimputed)                    
                            
               #normalizeBetweenArrays in limma
               #library(limma)
               normlog<-matrix(data = dataimputed, nrow =dim(dataimputed)[1], ncol = dim(dataimputed)[2], byrow = TRUE, dimnames = NULL)
               normlogbtwarray<-normalizeBetweenArrays(normlog,method='quantile')
               rownames(normlogbtwarray)=rownames(data)             
               colnames(normlogbtwarray)=colnames(data)

               return(normlogbtwarray)
}

#' Differentially expressed analysis
#' 
#' Find the top miRNAs and mRNAs that are differently expression between different conditions, e.g. cancer vs normal
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @param miR1 the miRNA dataset for condition 1, e.g. cancer
#' @param miR2 the miRNA dataset for condition 1, e.g. normal
#' @param mR1  the mRNA dataset for condition 1, e.g. cancer
#' @param mR2  the mRNA dataset for condition 2, e.g. normal
#' @param topkmiR the maximum number of miRNAs that we would like to extract, e.g. top 50 miRNAs.
#' @param topkmR the maximum number of mRNAs that we would like to extract, e.g. top 2000 mRNAs.
#' @param p.miR cutoff value for adjusted p-values when conducting differentially expressed analysis for miRNAs.
#' @param p.mR cutoff value for adjusted p-values when conducting differentially expressed analysis for mRNAs.
#' @return the dataset that includes differentially expressed miRNAs and mRNAs. columns are miRNAs and mRNAs  and rows are samples
#' @references
#' Smyth, G.K. (2005). Limma: linear models for microarray data. In Bioinformatics and computational biology solutions using R and Bioconductor (pp. 397-420). Springer New York.
#' @export
DiffExpAnalysis=function(miR1, miR2, mR1, mR2,topkmiR, topkmR, p.miR, p.mR){
	#miR1: miRNA dataset in condition 1, rows are miRNAs, samples are columns
	#miR2: condition2
	#mR1, mR2: genes in 2 conditions.
	# topk: the number of genes that we want to extract.
	
	#library(limma)
	miR1=read.csv(miR1, header=TRUE, sep=",")
	miRnames=miR1[,1]
	miR1=miR1[,-1]
	c1=ncol(miR1)
	
	miR2=read.csv(miR2, header=TRUE, sep=",")
	miR2=miR2[,-1]
	c2=ncol(miR2)
	
	miR=cbind(miR1, miR2)
	rownames(miR)=miRnames
	#miR=t(miR)
	
	mR1=read.csv(mR1, header=TRUE, sep=",")
	mRnames=mR1[,1]
	mR1=mR1[,-1]
	mR2=read.csv(mR2, header=TRUE, sep=",")
	mR2=mR2[,-1]
	mR=cbind(mR1, mR2)
	rownames(mR)=mRnames
	#mR=t(mR)
	Normal=NULL
  Cancer=NULL
	########## miR ############
	design=cbind(Normal=c(rep(1,c1), rep(0,c2)), Cancer=c(rep(0,c1), rep(1,c2)))
	miRfit=lmFit(miR, design)
	contrast.matrix=makeContrasts(NormalvCancer=Normal - Cancer, levels=design)
	miRfit2=contrasts.fit(miRfit, contrast.matrix)
	miRfit2=eBayes(miRfit2)
	miRresults=topTable(miRfit2, number= topkmiR, p.value=p.miR, sort.by="p", adjust.method="BH")
	write.csv(miRresults, file="DiffExpmiR.csv")
	########## mR ############
	mRfit=lmFit(mR, design)
	contrast.matrix=makeContrasts(NormalvCancer=Normal - Cancer, levels=design)
	mRfit2=contrasts.fit(mRfit, contrast.matrix)
	mRfit2=eBayes(mRfit2)
	mRresults=topTable(mRfit2, number= topkmR, p.value=p.mR, sort.by="p", adjust.method="BH")
	write.csv(mRresults, file="DiffExpmR.csv")
        miRSymbol=rownames(miRresults)
        mRSymbol=rownames(mRresults)
        miRExp=miR[which(miRnames %in% miRSymbol),]
        mRExp=mR[which(mRnames %in% mRSymbol),]
        
        miRExp=t(miRExp)
        mRExp=t(mRExp)
        
        DExp=cbind(miRExp,mRExp)
        
        return(DExp)
        	
}




#' miRNA target prediction with the Pearson correlation coefficient method
#' 
#' Calculate the Pearson correlation coefficient of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the Pearson correlation coefficients. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=Pearson(dataset, 1:3, 4:18) 
#' 
#' @references
#' Pearson, K. (1920) Notes on the history of correlation. Biometrika, 13, 25 - 45.
#' @export 

## 1. Pearson ##
Pearson=function(datacsv, cause, effect, targetbinding=NA){
data=Read(datacsv)
data=scale(data) #standardise the data
header<-readHeader(datacsv)
num_miRNA<-length(cause)
miR<-header[1:num_miRNA]
mR<-header[-(1:num_miRNA)]

miRNA=data[,cause]
mRNA=data[,effect]
corMat=cor(mRNA, miRNA, method="pearson")# miRNAs in columns and mRNAs in rows

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
corMat=corMat*edge
}

return(corMat)
}

#' miRNA target prediction with the Spearman correlation coefficient method
#' 
#' Calculate the Spearman correlation coefficient of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the Spearman correlation coefficients. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=Spearman(dataset, 1:3, 4:18) 
#' 
#' @references
#' Spearman, C. (1904) General intelligence, objectively determined and measured. Am. J. Psychol., 15, 201 - 92.
#' @export 
## 2. Spearman ##
Spearman=function(datacsv, cause, effect, targetbinding=NA){
data=Read(datacsv)
data=scale(data) #standardise the data
header<-readHeader(datacsv)
num_miRNA<-length(cause)
miR<-header[1:num_miRNA]
mR<-header[-(1:num_miRNA)]

miRNA=data[,cause]
mRNA=data[,effect]
corMat=cor(mRNA, miRNA, method="spearman")# miRNAs in columns and mRNAs in rows


if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
corMat=corMat*edge
}

return(corMat)
}

#' miRNA target prediction with the Kendall correlation coefficient method
#' 
#' Calculate the Kendall correlation coefficient of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the Kendall correlation coefficients. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=Kendall(dataset, 1:3, 4:18) 
#' @references
#' Kendall, M. (1938) A new measure of rank correlation. Biometrika, 30, 81 - 9.
#' @export 
## 3. Kendall ##
Kendall=function(datacsv, cause, effect, targetbinding=NA){
data=Read(datacsv)
data=scale(data) #standardise the data
header<-readHeader(datacsv)
num_miRNA<-length(cause)
miR<-header[1:num_miRNA]
mR<-header[-(1:num_miRNA)]

miRNA=data[,cause]
mRNA=data[,effect]
corMat=cor(mRNA, miRNA, method="kendall")# miRNAs in columns and mRNAs in rows


if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
corMat=corMat*edge
}

return(corMat)
}

#' miRNA target prediction with the Distance correlation  method
#' 
#' Calculate the Distance correlation  of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.
#' @importFrom energy dcov
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the Distance correlation values. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=Dcov(dataset, 1:3, 4:18) 
#' @references
#' Szekely, G., Rizzo, M. and Bakirov, N. (2007) Measuring and testing independence by correlation of distances. Ann. Stat., 35, 2769 - 94.
#' @export 
## 4. Dcov(Distance correlation) ##
Dcov <- function(datacsv, cause, effect, targetbinding=NA) {

       # library(energy)
        data=Read(datacsv)
        data=scale(data)
        header<-readHeader(datacsv)
        num_miRNA<-length(cause)
        miR<-header[1:num_miRNA]
        mR<-header[-(1:num_miRNA)]

        allnames=colnames(data)
	causenames=allnames[cause]
	effectnames=allnames[effect]
	Result<-matrix(nrow=length(effect), ncol=length(cause))
        for (i in cause){
                     for (j in effect){
                        Result[j-length(cause),i]<-dcov(data[,i],data[,j]) # Calculate Distance correlation.
                      }
         }
       colnames(Result)=causenames
       rownames(Result)=effectnames
       
       
if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
Result=Result*edge
}
 
return(Result)
}
#' miRNA target prediction with the Hoeffding correlation coefficient method
#' 
#' Calculate the Hoeffding correlation coefficient of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.
#' @importFrom Hmisc hoeffd 
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the Hoeffding correlation coefficients. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=Hoeffding(dataset, 1:3, 4:18) 
#' @references
#' Hoeffding, W. (1948) A non-parametric test of independence. Ann. Math. Stat., 19, 546 - 57.
#' @export 
## 5. Hoeffding(Hoeffding's D measure) ##
Hoeffding <- function(datacsv, cause, effect, targetbinding=NA) {

       #library(Hmisc)
        data=Read(datacsv)
        data=scale(data)
        header<-readHeader(datacsv)
        num_miRNA<-length(cause)
        miR<-header[1:num_miRNA]
        mR<-header[-(1:num_miRNA)]        

        allnames=colnames(data)
	causenames=allnames[cause]
	effectnames=allnames[effect]
	Result<-matrix(nrow=length(effect), ncol=length(cause))
        for (i in cause){
                     for (j in effect){
                        Result[j-length(cause),i]<-as.numeric(c(hoeffd(data[,i],data[,j])["D"], recursive = TRUE)[2]) # Calculate Hoeffding's D measure.
                      }
         }
       colnames(Result)=causenames
       rownames(Result)=effectnames

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
Result=Result*edge
}

return(Result)
}



#' @importFrom entropy mi.empirical
#' @importFrom gplots hist2d
      
EMI <- function(x, y) {
   # library(ghyp)
   # library(entropy)
    Mx <- matrix(x, length(x), 1)
    My <- matrix(y, length(y), 1)
    nbins <- ceiling(log(length(Mx[,1]),2)) + 1
    a <- hist2d(cbind(Mx,My), nbins=nbins,show=FALSE)
    result <- mi.empirical(a$counts)
    return(result)
}

#' miRNA target prediction with  mutual information method
#' 
#' Calculate the mutual information of each pair of miRNA-mRNA,and return a matrix of mutual information values with columns are miRNAs and rows are mRNAs.
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the mutual information values. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=MI(dataset, 1:3, 4:18) 
#' @references
#' Moon, Y.I., Balaji, R., and Lall, U. (1995) Estimation of mutual information using kernel density estimators. Phys. Rev. E, 52, 2318 - 21.
#' @export 

MI <- function(datacsv, cause, effect, targetbinding=NA) {
        data=Read(datacsv)
        data=scale(data)
        header<-readHeader(datacsv)
        num_miRNA<-length(cause)
        miR<-header[1:num_miRNA]
        mR<-header[-(1:num_miRNA)]

        allnames=colnames(data)
	causenames=allnames[cause]
	effectnames=allnames[effect]
	Result<-matrix(nrow=length(effect), ncol=length(cause))
        for (i in cause){
                     for (j in effect){
                        Result[j-length(cause),i]<-EMI(data[,i],data[,j]) # Calculate Mutual Information.
                      }
         }
       colnames(Result)=causenames
       rownames(Result)=effectnames       

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
Result=Result*edge
}

return(Result)
}

#' miRNA target prediction with the IDA method
#' 
#' Calculate the causal effect of each pair of miRNA-mRNA,and return a matrix of causal effects with columns are miRNAs and rows are mRNAs.
#' @importFrom pcalg pc idaFast gaussCItest
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param pcmethod choose different versons of the PC algorithm, including "original" (default)
#' "stable", and "stable.fast"
#' @param alpha significance level for the conditional independence test, e.g. 0.05.
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the causal effects. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=IDA(dataset, 1:3, 4:18) 
#' @references
#' 1. Le, T.D., Liu, L., Tsykin, A., Goodall, G.J., Liu, B., Sun, B.Y. and Li, J. (2013) Inferring microRNA-mRNA causal regulatory relationships from expression data. Bioinformatics, 29, 765-71.
#' 
#' 2. Zhang, J., Le, T.D., Liu, L., Liu, B., He, J., Goodall, G.J. and Li, J. (2014) Identifying direct miRNA-mRNA causal regulatory relationships in heterogeneous data. J. Biomed. Inform., 52, 438-47.
#' 
#' 3. Maathuis, H.M., Colombo, D., Kalisch, M. and Buhlmann, P. (2010) Predicting causal effects in large-scale systems from observational data. Nat. Methods, 7, 247-249.
#' 
#' 4. Maathuis, H.M., Kalisch, M. and Buhlmann, P. (2009) Estimating high-dimensional intervention effects from observational data. Ann. Stat., 37, 3133-3164.
#' @export 
## 7. IDA ##
IDA=function(datacsv, cause, effect, pcmethod="original", alpha=0.05, targetbinding=NA){
     #   library(pcalg)

	if(is.character(datacsv)){
		data=Read(datacsv)
		#data=data[,-1] # if the dataset have the sample names column, otherwise comment this out.
	} else {
		data=datacsv #Assume there is no samplenames column and this is a data.frame.
	}							
	data=scale(data) #standardise the data
	header<-readHeader(datacsv)
        num_miRNA<-length(cause)
        miR<-header[1:num_miRNA]
        mR<-header[-(1:num_miRNA)]

        #print(data[1:5,])
	allnames=colnames(data)
	causenames=allnames[cause]
	effectnames=allnames[effect]
		
	multiset=character(0)
	result=matrix(nrow=length(effect), ncol=length(cause))
	suffStat=list(C=cor(data), n=nrow(data))
	indepTest=gaussCItest
	
	pcFit <- pc(suffStat, indepTest, p=ncol(data), alpha=alpha, skel.method=pcmethod)
	for (l in cause){
				
			#Inferring causal effects
			caef<-idaFast(l,effect,cov(data), pcFit@graph )
		
			#min of absolute values.
			caef1<-matrix(nrow=length(effect),ncol=1)
			for (k in 1:length(effect)){
				caefabs<-abs(caef)
				index<-which(caefabs==min(caefabs[k,]), arr.ind=TRUE)
				pos<-index[1,2]
				caef1[k,]<-caef[k,pos]
			}
			result[,l]<-caef1
	}
	colnames(result)=causenames
	rownames(result)=effectnames

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
result=result*edge
}

return(result)	
}



## 9. RDC(Randomized Dependence Coefficient) ##
RDCParameter <- function(x,y,k=20,s=1/6,f=sin) {
  x <- cbind(apply(as.matrix(x),2,function(u)rank(u)/length(u)),1)
  y <- cbind(apply(as.matrix(y),2,function(u)rank(u)/length(u)),1)
  x <- s/ncol(x)*x%*%matrix(rnorm(ncol(x)*k),ncol(x))
  y <- s/ncol(y)*y%*%matrix(rnorm(ncol(y)*k),ncol(y))
  cancor(cbind(f(x),1),cbind(f(y),1))$cor[1]
}

#' miRNA target prediction with the Randomized Dependence Coefficient method
#' 
#' Calculate the Randomized Dependence coefficient of each pair of miRNA-mRNA,and return a matrix of  coefficients with columns are miRNAs and rows are mRNAs.
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the  correlation coefficients. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=RDC(dataset, 1:3, 4:18) 
#' @export 

RDC <- function(datacsv,cause, effect, targetbinding=NA) {
        data=Read(datacsv)
        data=scale(data) #standardise the data
        header<-readHeader(datacsv)
        num_miRNA<-length(cause)
        miR<-header[1:num_miRNA]
        mR<-header[-(1:num_miRNA)]

        allnames=colnames(data)
	causenames=allnames[cause]
	effectnames=allnames[effect]
        
	Result<-matrix(nrow=length(effect), ncol=length(cause))
        for (i in cause){
                     for (j in effect){
                        Result[j-length(cause),i]<-RDCParameter(data[,i],data[,j],k=20,s=1/6,f=sin) # Calculate Randomized Dependence Coefficient.
                      }
         }

        colnames(Result)=causenames
	rownames(Result)=effectnames

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
Result=Result*edge
}
        
return(Result)
}
#' miRNA target prediction with the Lasso method
#' 
#' Calculate the Lasso regression coefficient of each pair of miRNA-mRNA, and return a matrix of coefficients with columns are miRNAs and rows are mRNAs.
#' @importFrom glmnet glmnet
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the Lasso regression coefficients. Columns are miRNAs, rows are mRNAs.
#' @examples
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=Lasso(dataset, 1:3, 4:18) 
#' @references
#' 1. Le, T.D., Zhang, J., Liu, L., and Li, J. (2015) Ensemble Methods for miRNA Target Prediction from Expression Data, submitted.
#' 
#' 2. Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. J. R. Stat. Soc. Series B Stat. Methodol., 267-288.
#' @export 
## 10. Lasso ##
Lasso=function(datacsv, cause, effect, targetbinding=NA){
 # library(glmnet)
	if(is.character(datacsv)){
		data=Read(datacsv)
		#data=data[,-1] # if the dataset have the sample names column, otherwise comment this out.
	} else {
		data=datacsv #Assume there is no samplenames column and this is a data.frame.
	}							#To allow both .csv data input or a matrix in R. This will help the IDAbootstrap, as IDA can run on sampling matrices.
	data=scale(data) #standardise the data
        header<-readHeader(datacsv)
        num_miRNA<-length(cause)
        miR<-header[1:num_miRNA]
        mR<-header[-(1:num_miRNA)]
  data=as.matrix(data)
	#print(data[1:5,])
	allnames=colnames(data)
	causenames=allnames[cause]
	effectnames=allnames[effect]
	res=matrix(nrow=length(effectnames), ncol=length(causenames), 0)
	colnames(res)=causenames
	rownames(res)=effectnames
		
	for(gene in effectnames) {
		aGene <- glmnet(data[,cause],data[,gene],alpha=1)$beta #return the the effects of all miRNAs on the gene
    aGene=as.matrix(aGene)
		aGene=rowMeans(aGene) # take the means of all values output from lasso    
		res[gene,]=aGene # assign the effects of all miRNAs on the gene to the result. 
	}

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
res=res*edge
}

return(res)
}

#' miRNA target prediction with the Elastic-net regression coefficient method
#' 
#' Calculate the Elastic-net regression coefficient of each pair of miRNA-mRNA,and return a matrix of correlation coefficients with columns are miRNAs and rows are mRNAs.
#' @importFrom glmnet glmnet
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the Elastic-net regression coefficients. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=Elastic(dataset, 1:3, 4:18) 
#' @references
#' 1. Le, T.D., Zhang, J., Liu, L., and Li, J. (2015) Ensemble Methods for miRNA Target Prediction from Expression Data, under review.
#' 
#' 2. Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. J. R. Stat. Soc. Series B Stat. Methodol., 67, 301-320.
#' @export 
## 11. Elastic-net ##
Elastic=function(datacsv, cause, effect, targetbinding=NA){
 # library(glmnet)
	if(is.character(datacsv)){
		data=Read(datacsv)
		#data=data[,-1] # if the dataset have the sample names column, otherwise comment this out.
	} else {
		data=datacsv #Assume there is no samplenames column and this is a data.frame.
	}							#To allow both .csv data input or a matrix in R. This will help the IDAbootstrap, as IDA can run on sampling matrices.
	data=scale(data) #standardise the data
        header<-readHeader(datacsv)
        num_miRNA<-length(cause)
        miR<-header[1:num_miRNA]
        mR<-header[-(1:num_miRNA)]
	data=as.matrix(data)
	#print(data[1:5,])
	allnames=colnames(data)
	causenames=allnames[cause]
	effectnames=allnames[effect]
	res=matrix(nrow=length(effectnames), ncol=length(causenames), 0)
	colnames(res)=causenames
	rownames(res)=effectnames
		
	for(gene in effectnames) {
		aGene <- glmnet(data[,cause],data[,gene],alpha=0.5)$beta #return the the effects of all miRNAs on the gene
		aGene=as.matrix(aGene)
    aGene=rowMeans(aGene) # take the means of all values output from lasso     
		res[gene,]=aGene # assign the effects of all miRNAs on the gene to the result. 
	}

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
res=res*edge
}

return(res)
}

#' miRNA target prediction with the Z-score method
#' 
#' Calculate the Z-score value of each pair of miRNA-mRNA, and return a matrix of values with columns are miRNAs and rows are mRNAs.
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the Z-score values. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=Zscore(dataset, 1:3, 4:18) 
#' @references
#' Prill, R.J., Marbach, D., Saez-Rodriguez, J., Sorger, P.K., Alexopoulos, L.G., Xue, X., Clarke, N.D., Altan-Bonnet, G. and Stolovitzky, G. (2010) Towards a rigorous assessment of systems biology models: the DREAM3 challenges. PLoS One, 5, e9202.
#' @export 
## 12. Z-score ##
Zscore=function(datacsv, cause, effect, targetbinding=NA){
	dt=read.csv(datacsv, header=TRUE, sep=",")
	dt=scale(dt)
        header<-readHeader(datacsv)
        num_miRNA<-length(cause)
        miR<-header[1:num_miRNA]
        mR<-header[-(1:num_miRNA)]

	effectnames=colnames(dt)[effect]
	causenames=colnames(dt)[cause]
	res=matrix(nrow=length(effectnames), ncol=length(causenames), 0)
	colnames(res)=causenames
	rownames(res)=effectnames
	causes=dt[,cause]
	effects=dt[,effect]
	#zAB=(xBminA-meanB)/sdB
	for (i in 1:length(cause)){
		for (j in 1: length(effect)){
			indexminA=which(causes[,i]==min(causes[,i]))
			xBminA=effects[indexminA,j]
			xBminA=median(xBminA)
			#sdB=sd(effects[,j]) #if we standardise the sdB=1
			#meanB=mean(effects[,j]) if we standardise the mean is 0
			#zij=abs(xBminA-meanB)/sdB
			zij=abs(xBminA)
			
			res[j,i]=zij
		}
	}

if(is.na(targetbinding)==FALSE){
#query knowledge matrix from file
edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
edgeTargetScan<-edgeTargetScan!=0;
edge=edgeTargetScan[effect,cause]
res=res*edge
}

return(res)
}

#' miRNA target prediction with the ProMISe method
#' 
#' Calculate the ProMISe score of each pair of miRNA-mRNA, and return a matrix of values with columns are miRNAs and rows are mRNAs.
#' @importFrom Roleswitch roleswitch
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param targetbinding the putative target, e.g. "TargetScan.csv". If targetbinding is not specified, only expression data is used.
#' If targetbinding is specified, the prediction results using expression data with be intersected with the interactions in the target binding file.
#' @return A  matrix that includes the ProMISe scores. Columns are miRNAs, rows are mRNAs.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' results=ProMISe(dataset, 1:3, 4:18) 
#' @references
#' Li, Y., Liang, C., Wong, K.C., Jin, K., and Zhang, Z. (2014). Inferring probabilistic miRNA - mRNA interaction signatures in cancers: a role-switch approach. Nucleic Acids Res., 42, e76-e76.
#' @export 
## 13. ProMISe ##
ProMISe=function(datacsv, cause, effect, targetbinding=NA){
  # library("Roleswitch")
  dt<-Read(datacsv)
  dd<-colMeans(dt)
  stdData<-as.matrix(dd)
  header<-readHeader(datacsv)
  num_miRNA<-length(cause)
  miR<-header[1:num_miRNA]
  mR<-header[-(1:num_miRNA)]
  
  x.o<-matrix(stdData[effect,],dimnames=list(c(1:length(effect)),"mRNA"))
  z.o<-matrix(stdData[cause,],dimnames=list(c(1:length(cause)),"miRNA"))
  c<-matrix(1,length(effect),length(cause)) #Generate ones matrix
  rownames(c)<-c(1:length(effect))
  colnames(c)<-c(1:length(cause))
  
  rMatrix <- roleswitch(x.o,z.o,c)$p.xz # Calculate ProMISe probabilistic
  rownames(rMatrix) <- colnames(dt)[effect]
  colnames(rMatrix) <- colnames(dt)[cause]
  
  if(is.na(targetbinding)==FALSE){
    #query knowledge matrix from file
    edgeTargetScan<-queryTargetFile(miR,mR,targetbinding); 
    edgeTargetScan<-edgeTargetScan+t(edgeTargetScan);
    edgeTargetScan<-edgeTargetScan!=0;
    edge=edgeTargetScan[effect,cause]
    rMatrix=rMatrix*edge
  }
  
  return(rMatrix)
}


#' Extract top k miRNA-mRNA interactions 
#' 
#' Rank the miRNA-mRNA interactions based on absolute values of the correlations/scores/causal effects, and return
#' the topk interactions.
#' @param cormat the correlation matrix that need to be extracted with columns are miRNAs and rows are mRNAs
#' @param topk the number of interactions that need to be extracted.
#' @return topk interactions 
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' EMTresults=Pearson(dataset, 1:3, 4:18)
#' top10=Extopk(EMTresults, 10)
#' @export
Extopk <- function(cormat,topk){
          mir<-ncol(cormat)
          gene<-nrow(cormat)
          mirname<-colnames(cormat)
          mirname = gsub("\\.", "-", mirname)  # replace "." with "-" for miRNA
          for(index in 1:length(mirname)){
          if(substr(mirname, nchar(mirname),nchar(mirname))[index]=="-") {
	          substring(mirname[index], nchar(mirname[index]), nchar(mirname[index]))="*" 
                  }
          }
          genename<-rownames(cormat)
          Result<-matrix(-1,ncol=4,nrow=mir*gene)
          n<-0
          for(i in 1:mir){
                   for(j in 1:gene){
                       n<-(n+1)
                       Result[n, 1] <- mirname[i]
                       Result[n, 2] <- genename[j]
                       Result[n, 3] <- cormat[j,i]
                       Result[n, 4] <- abs(as.numeric(Result[n, 3])) # Calculate absolute value
                }
           } 
           TrResult <- Result[sort.list(as.numeric(Result[,4]),decreasing=TRUE),] #Rank the whole interations by decreasing the absolute value
           
           return(TrResult[1:topk,]) # Extract top k miRNA-mRNA interactions
        	
}      


#' Ensemble method for miRNA target prediction using Borda count election
#' 
#' Use the Borda count election method to integrate the rankings from different miRNA target prediction methods
#' @param listCEmatrices a list of matrices that include the correlation coefficients/causal effects/scores resulting from different target prediction methods
#' @return a matrix of ranking scores (averaging all the rankings from different methods). Columns are miRNAs and rows are mRNAs
#' @examples
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' ps=Pearson(dataset, cause=1:3, effect=4:18)
#' ida=IDA(dataset, cause=1:3, effect=4:18)
#' borda=Borda(list(ps, ida))
#' @references
#' 1. Le, T.D., Zhang, J., Liu, L., and Li, J. (2015) Ensemble Methods for miRNA Target Prediction from Expression Data, Plos ONE.
#' 
#' 2. Marbach, D., Costello, J.C., Kuffner, R., Vega, N.M., Prill, R.J., Camacho, D.M., Allison, K.R. and DREAM5 Consortium (2012). Wisdom of crowds for robust gene network inference. Nat. Methods, 9, 796-804.
#' @export

Borda=function(listCEmatrices){
noMethods=length(listCEmatrices)
effectnames=rownames(listCEmatrices[[1]])
causenames=colnames(listCEmatrices[[1]])
res=matrix(nrow=length(effectnames), ncol=length(causenames), 0)
colnames(res)=causenames
rownames(res)=effectnames

for(i in 1:noMethods){ # for each CEmatrix from a method
	for (j in 1:length(causenames)){ #for each miRNA
		cormat=listCEmatrices[[i]][,j] #extract the corelation values
		cormat=cormat[order(-abs(cormat))] # sort based on absolute values
		for(k in 1:length(cormat)){cormat[k]=k} # change the correlation values by the ranking
		rn=names(cormat) # take the genenames
		
		for(gene in rn){ #for each gene in the sorted matrix
			res[gene, j]=res[gene, j]+cormat[gene] #add up the current rankings of the possition (gene, miRNA)
		}
	}
}
res=res/noMethods #take the average rankings
res=length(effectnames)/res #revert the rankings as bRank will sort the rankings as if correlation.
}

#' Ensemble method for miRNA target prediction using Borda count election with topk targets
#' 
#' Use the Borda count election method to integrate the rankings from different miRNA target prediction methods, but only topk targets of each miRNA are included
#' in the calculation. The targets outside the topk will be assigned a large and fixed rank, e.g. number of genes in the dataset.
#' @param listCEmatrices a list of matrices that include the correlation/causal effects/scores resulting from a target prediction method
#' @param topk number of targets of a miRNA to be included in the calculation (Borda count election)
#' @return a matrix of ranking scores (averaging all the rankings from different methods). Columns are miRNAs and rows are mRNAs
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' ps=Pearson(dataset, cause=1:3, effect=4:18)
#' ida=IDA(dataset, cause=1:3, effect=4:18)
#' borda=BordaTopk(list(ps, ida), topk=10)
#' @references
#' Le, T.D., Zhang, J., Liu, L., and Li, J. (2015) Ensemble Methods for miRNA Target Prediction from Expression Data, Plos ONE.
#' @export
BordaTopk=function(listCEmatrices, topk){
noMethods=length(listCEmatrices)
effectnames=rownames(listCEmatrices[[1]])
causenames=colnames(listCEmatrices[[1]])
res=matrix(nrow=length(effectnames), ncol=length(causenames), 0)
colnames(res)=causenames
rownames(res)=effectnames
noGenes=nrow(listCEmatrices[[1]])
for(i in 1:noMethods){ # for each CEmatrix from a method
	for (j in 1:length(causenames)){ #for each miRNA
		cormat=listCEmatrices[[i]][,j] #extract the corelation values
		cormat=cormat[order(-abs(cormat))] # sort based on absolute values
		for(k in 1:length(cormat)){cormat[k]=k} # change the correlation values by the ranking
		rn=names(cormat) # take the genenames
		
		for(gene in rn){ #for each gene in the sorted matrix
			if(cormat[gene]>topk) cormat[gene]=noGenes
			res[gene, j]=res[gene, j]+cormat[gene] #add up the current rankings of the possition (gene, miRNA)
		}
	}
}
res=res/noMethods #take the average rankings
res=length(effectnames)/res #revert the rankings as bRank will sort the rankings as if correlation.
}



ReOrder=function(topkList){
   
   ####### preprocessing the list of topk results #######
   topkList = as.matrix(topkList, ncol=3);
   if(nrow(topkList)<1 || length(topkList[,3]>0)==0 || length(topkList[,3]<0)==0)   return(data.frame(topkList));
   
   ####### sort the values in increasing order ########
   val = as.numeric(topkList[,3]);  
   val = sort(val, index.return=TRUE, decreasing=FALSE)

   ####### obtain the negative values ########
   ni = val$ix[which(val$x<0)] # the index of negative values in Val;
   pi = val$ix[length(val$x):(length(ni)+1)] # the index of positive values in Val

   return(data.frame(topkList[c(ni, pi), ]))
}


#' Extract topk predicted targets of a miRNA
#' Rank all the targets of a miRNA and extract the topk targets
#' @param CEmatrix the matrix of correlation/causal effect/score results with columns are miRNAs and rows are mRNAs
#' @param causeIndex the column index of the miRNA that we would like to extract
#' @param topk the number of targets being extracted
#' @param downreg if TRUE the negative correlation/causal effect/score will be on the top of the ranking. This is to
#' favour the negative regulations. 
#' @return a matrix with 3 columns, where the first column contains the miRNA, the second column contains the mRNAs and the last column contains the correlations/causal effects/scores
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' ps=Pearson(dataset, cause=1:3, effect=4:18)
#' miR200aTop10 = bRank(ps, 3, 10, TRUE)
#' @export

bRank=function(CEmatrix, causeIndex, topk, downreg=TRUE){

miRNAnames=colnames(CEmatrix)
mRNAnames=rownames(CEmatrix)

c1=rep(miRNAnames[causeIndex], length(mRNAnames))
corMat=CEmatrix[,causeIndex]
result=cbind(c1,mRNAnames, corMat)
rownames(result)=NULL
colnames(result)=c("miRNA", "mRNA", "Correlation")
#remove quotes
result=data.frame(result)
result[,3]=as.numeric(as.character(result[,3]))

#cat("numberof genes in result:", nrow(result), "\n")
result=result[order(-(abs(result[,3]))),]
result=result[!duplicated(result[,2]),]
#cat("number of genes in the ranking list after removing dulicates", nrow(result), "\n")

result=result[1:topk,]
result[, 1] = gsub("\\.", "-", result[, 1])  # replace "." with "-" for miRNA
if(substr(result[ ,1], nchar(result[,1]), nchar(result[,1]))[1]=="-") {
	substring(result[,1], nchar(result[,1]), nchar(result[,1]))="*" 
}
# if the last character is - then change to *, e.g. hsa-miR-7-1*
#result[, 2] = gsub("\\.", "-", result[, 2])  # replace "." with "-" for mRNA
if(downreg) result=ReOrder(result)
#print(result[1:5,])
return(result)
}


############# Validation function ###### 
#' Validate the targets of a miRNA
#' 
#' Given the predicted target of a miRNA, the function returns a list of targets that are experimentally confirmed
#' based on the provided ground truth. Users can provide their own ground truth or use the built-in ground truth 
#' which is the union of Tarbase, miRTarbase, miRecords, and miRWalk. 
#' @param topkList a matrix with 3 columns. The first column is the miRNA name, the second contains the target mRNAs, and the third contains the correlation values/ causal effects/ scores
#' @param datacsv the ground truth for the validation. The ground truth is a matrix with 2 columns, where the first column
#' is the miRNA and the second is the mRNA.
#' @return a matrix in the same format of the input matrix put only contains the confirmed interactions.
#' @examples 
#' dataset=system.file("extdata", "ToyEMT.csv", package="miRLAB")
#' ps=Pearson(dataset, cause=1:3, effect=4:18)
#' miR200aTop10=bRank(ps, 3, 10, TRUE)
#' groundtruth=system.file("extdata", "Toygroundtruth.csv", package="miRLAB")
#' miR200aTop10Confirmed = Validation(miR200aTop10, groundtruth)
#' @export
Validation=function(topkList, datacsv){
	
####### preprocessing the list of topk results #######
	topkList = as.matrix(topkList, ncol=3);
	if(nrow(topkList)<1) stop("No data in the input")
    
####### read the validation data from file ########
        dt=read.csv(datacsv, header=TRUE, sep=",")
	
####### get the validation lists from the data ######
	        
        dt = paste(dt[, 1], dt[, 2], sep=" ");
        tmp= paste(topkList[, 1], topkList[, 2], sep=" ");	  	
	result=topkList[which(tmp %in% dt), ] 
        
	
	if(is.matrix(result)){
		result=data.frame(result)
		result[,3]=as.numeric(as.character(result[,3]))
		NoConfimred=nrow(result)
	} else {
		result=as.matrix(result)
		result=t(result)
		result=data.frame(result)
		result[,3]=as.numeric(as.character(result[,3]))
		NoConfimred=nrow(result)
	}
	
	return(list(result, NoConfimred))
}

#' Validate the targets of a miRNA using transfection data
#' 
#' Given the predicted target of a miRNA, the function returns a list of targets that are  confirmed
#' based on the curated transfection data. Users need to download the file logFC.imputed.rda from nugget.unisa.edu.au/Thuc/miRLAB/ and place it in the working directory (this file is obtained from the TargetScoreData package)
#' @param topkList a matrix with 3 columns. The first column is the miRNA name, the second contains the target mRNAs, and the third contains the correlation values/ causal effects/ scores
#' @param LFC the log fold change threshold. The targets that have the absolute value of log fold change greater than the LFC
#' will be regarded as the confirmed targets.
#' @return a matrix in the same format of the input matrix put only contains the confirmed interactions.
#' @examples 
#dataset=system.file("extdata", "EMT35.csv", package="miRLAB")
#' print("ps=Pearson(dataset, cause=1:35, effect=36:1189)")
#' print("miR200aTop100=bRank(ps, 11, 100, TRUE)")
#' print("miR200aTop100Confirmed = ValidationT(miR200aTop100, 1.0)")
#' @export
#' @references
#' 1. Le, T.D., Zhang, J., Liu, L., and Li, J. (2015) Ensemble Methods for miRNA Target Prediction from Expression Data, under review.
#' 
#' 2. Li Y, Goldenberg A, Wong K and Zhang Z (2014). A probabilistic approach to explore human microRNA targetome using microRNA-overexpression data and sequence information. Bioinformatics, 30(5), pp. 621-628. http://dx.doi.org/10.1093/bioinformatics/btt599.


ValidationT=function(topkList, LFC){
 #get all transfection data from TargetScoreData package and present in the format of Tarbase.
 #the rows with LFC <= LFC will be discarded
 dt=system.file("extdata", "Transfection_data_summary.csv", package="miRLAB")
 Transfection=read.csv(dt, header=TRUE, sep=",")
 logFC.imputed=c()
 if(!file.exists("logFC.imputed.rda"))
   stop("Please download the transfection data from nugget.unisa.edu.au/Thuc/miRLAB/logFC.imputed.rda for validation. 
       If you would like to use experimentally confirmed data only, please use the Validation function")
 load("./logFC.imputed.rda")
 td=logFC.imputed
 
 topkList = as.matrix(topkList, ncol=3);
 stopifnot(nrow(topkList)>0)
 miN = unique(topkList[, 1]) # get the unique names of miRNA
 mind= NULL;
 for(i in 1:length(miN))
	mind = c(mind, which(Transfection[,1]==miN[i]))  #get the index of miRNA from the validation data   
 #print(mind)
 if(length(mind)==0){
	#cat("This miRNA is not in the tranfection database", "\n")
	return(list(NULL, 0))
}	
 rn=rownames(td) # take the genenames
 td = td[,mind]      # get a small validation data from the original one, a new list
 
 if(length(mind)>1) td=rowMeans(td)	 # take the average of LFCs
 		 
#decorate the table
 groundtruth=cbind(rep(miN, length(rn)), rn, td)
 rownames(groundtruth)=NULL
 colnames(groundtruth)=c("miRNA", "mRNA", "LFC")
 groundtruth=data.frame(groundtruth)
 groundtruth[,3]=as.numeric(as.character(groundtruth[,3]))
 groundtruth=groundtruth[abs(groundtruth[,3])>LFC,]       #LFC<-1 using groundtruth=groundtruth[groundtruth[,3]<-LFC,]
 groundtruth=groundtruth[,1:2]
 
 groundtruth = paste(groundtruth[, 1], groundtruth[, 2], sep=" ");
 tmp= paste(topkList[, 1], topkList[, 2], sep=" ");
 result=topkList[which(tmp %in% groundtruth), ]
 
 #Decorate the result table, if only one row it is not a matrix. If it is a matrix we just remove the "" in data
 if(is.matrix(result)){
	result=data.frame(result)
	result[,3]=as.numeric(as.character(result[,3]))
	NoConfimred=nrow(result)
  } else {
		result=as.matrix(result)
		result=t(result)
		result=data.frame(result)
		result[,3]=as.numeric(as.character(result[,3]))
		NoConfimred=nrow(result)
  }
  #cat("Number of confirmed genes for  ",miN, "is: ", NoConfimred, "\n")
	
  return(list(result, NoConfimred))

 
}

#' Validate the targets of all miRNA using both experimentally confirmed and transfection data
#' 
#' Given the predicted target of all miRNA, the function returns a list of targets of each miRNA that are  confirmed
#' based on the experimentally validated interactions or curated transfection data. Users need to download the file logFC.imputed.rda from nugget.unisa.edu.au/Thuc/miRLAB/ and place it in the working directory (this file is obtained from the TargetScoreData package)
#' @param CEmatrix the matrix of correlation/causal effects/scores with columns are miRNAs and rows are mRNAs
#' @param topk the number of targets of each miRNA that are being validated. 
#' @param groundtruth the csv file containing the ground truth.
#' @param LFC the log fold change threshold for the transfection data. The targets that have the absolute value of log fold change greater than the LFC
#' will be regarded as the confirmed targets.
#' @param downreg if TRUE the negative correlation/causal effect/score values will be ranked on the top of the ranking. This is to favour the down regulations.
#' @return a list of matrices that contains the confirmed interactions by both provided ground truth and built-in transfection data.
#' @examples 
#' print("ps=Pearson(dataset, cause=1:3, effect=4:18)")
#' print("results=ValidateAll(ps, 10, groundtruth, LFC=0.5, downreg=TRUE)")
#' @export
# \dontrun{
# dataset=system.file("extdata", "ToyEMT35.csv", package="miRLAB")
# ps=Pearson(dataset, cause=1:3, effect=4:18)
# groundtruth=system.file("extdata", "groundtruth.csv", package="miRLAB")
# results=ValidateAll(ps, 10, groundtruth, LFC=0.5, downreg=TRUE)
# }
ValidateAll=function(CEmatrix, topk, groundtruth, LFC, downreg=TRUE){
#CEmatrix: the results from a computational method in a matrix format. columns are miRNAs. Rows are genes.
#causes: the column indices of the causes in the dataset or in the CEMatrix.
#Top k gene of each miRNA we would like to extract for validation.
#Groundtruth is the experimentally validated database in .csv format.
#LFC: log2 fold-change threshold for identifying the groundtruth using miRNA transfection data.
if(!file.exists("logFC.imputed.rda"))
  stop("Please download the transfection data from nugget.unisa.edu.au/Thuc/miRLAB/logFC.imputed.rda for validation. 
       If you would like to use experimentally confirmed data only, please use the Validation function")
 causes=1:ncol(CEmatrix)
 NoExpConfirmed=c()
 NoTransConfirmed=c()
 names=c('miRNA','mRNA','Correlation')
 K1=c() #ground truth experiment
 ResultK1=matrix(,ncol=3)
 colnames(ResultK1)=names
 K2=c() #ground truth transfection
 ResultK2=matrix(,ncol=3)
 colnames(ResultK2)=names
 pE=c() #pvalue for experiment
 pT=c() #pvalue for transfection
 S=nrow(CEmatrix)
	for (i in causes){
		miRtopk=bRank(CEmatrix, i, topk, downreg)
		miRall=bRank(CEmatrix, i, S)
		#cat("Validate the prediction using experimentally confirmed database ", "\n")
		temp1=Validation(miRtopk, groundtruth)
		temp2=Validation(miRall, groundtruth)
		pvalE=1-phyper((temp1[[2]]-1), temp2[[2]], (S-temp2[[2]]), topk)
		pE=c(pE, pvalE)
		
		cat("EXPERIMENT:  ", colnames(CEmatrix)[i], ": ", temp1[[2]], "with ", temp2[[2]], " in the groundtruth. pvalue: ",pvalE, "\n")
		NoExpConfirmed=c(NoExpConfirmed,temp1[[2]]) # sum will be x in hypergeometric test.
		K1=c(K1,temp2[[2]])
                if(temp1[[2]]>0)  ResultK1=rbind(ResultK1,temp1[[1]])
                     
 		
		#cat("Validate the prediction using transfection data ", "\n")
		temp3=ValidationT(miRtopk, LFC)
		temp4=ValidationT(miRall, LFC)
		pvalT=1-phyper((temp3[[2]]-1), temp4[[2]], (S-temp4[[2]]), topk)
		pT=c(pT, pvalT)
		cat("TRANSFECTION:  ", colnames(CEmatrix)[i], ": ", temp3[[2]], "with ", temp4[[2]], " in the groundtruth. pvalue: ", pvalT, "\n")
		NoTransConfirmed=c(NoTransConfirmed,temp3[[2]]) # sum will be x in the hypergeometric test.
		K2=c(K2,temp4[[2]]) 
                
		if(temp3[[2]]>0) ResultK2=rbind(ResultK2,temp3[[1]])
		
	}
	pvalueE=1-phyper((sum(NoExpConfirmed)-1), sum(K1), (S*ncol(CEmatrix)-sum(K1)), topk*ncol(CEmatrix))
	pvalueT=1-phyper((sum(NoTransConfirmed)-1), sum(K2), (S*ncol(CEmatrix)-sum(K2)), topk*ncol(CEmatrix))
	
	cat("number of confirms by experiments: ", sum(NoExpConfirmed), ", pvalue: ", pvalueE, ". By transfection ", sum(NoTransConfirmed), ", pvalue: ", pvalueT, "\n")
	cat("groundtruths in experiments: ", sum(K1), ", and in transfection ", sum(K2), "\n")
	result=list(ResultK1[2:nrow(ResultK1),], cbind(K1, NoExpConfirmed, pE), ResultK2[2:nrow(ResultK2),], cbind(K2, NoTransConfirmed, pT))
}


## Read external results ##
#' Read results from other methods
#' 
#' Read the results predicted by external methods (methods that are not in this package and may not be implemented in R). Consequently, we can compare the results
#' predicted by the external methods and results predicted by the methods in the miRLAB package.
#' @param datacsv the input dataset in csv format
#' @param cause the column range that specifies the causes (miRNAs), e.g. 1:35
#' @param effect the column range that specifies the effects (mRNAs), e.g. 36:2000
#' @param ExtCEcsv score matrix predicted by an external matrix with columns are miRNAs and rows are mRNAs.
#' @return a matrix of scores predicted by an external matrix and ready for further validation and comparison tasks.
#' @examples 
#' print("GenemiR=ReadExtResult(dataset, cause=1:3, effect=4:18, 'genemirresults.csv')")
#' @export
ReadExtResult=function(datacsv, cause, effect,  ExtCEcsv){
	dataset=read.csv(datacsv, header=TRUE, sep=",")
	genenames=colnames(dataset)
	causenames=genenames[cause]
	effectnames=genenames[effect]
	Extmatrix=read.csv(ExtCEcsv, header=FALSE, sep=",")
	colnames(Extmatrix)=causenames
	rownames(Extmatrix)=effectnames
	return(Extmatrix)
	
}

#delete experiment function

## Compare the validation results of 13 built-in methods ##
filterAndCompare=function(allresults, noVal){
	#allresults: the results from all methods generated from experiment function. This is a list
	#noVal: number of confirmed targets in each method (threshold) to filter. Records (miRNA) with less than this will be removed
	ExpResult=allresults[[1]]
	TransResult=allresults[[2]]
	temp1=apply(ExpResult, 1, function(x) all(x[c(2,4,6,8,10,12,14,16,18,20,22,24,26)]>(noVal-1)))
	ExpResult=ExpResult[temp1,]
	temp2=apply(TransResult, 1, function(x) all(x[c(2,4,6,8,10,12,14,16,18,20,22,24,26)]>(noVal-1)))
	TransResult=TransResult[temp2,]
	###If else to incase only one record. In that case the result is not a matrix
	if(is.matrix(ExpResult)){
		ExpResult=ExpResult[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
		colnames(ExpResult)=c( "Pearson", "Spearman", "Kendall", " Dcov", "Hoeffding", "MI", "IDA", "MIC", "RDC", "Lasso", "Elastic", "Z-score", "ProMISe")
	} else {
		tt=as.matrix(ExpResult)
		ExpResult=t(tt)
		ExpResult=ExpResult[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
		names(ExpResult)=c( "Pearson", "Spearman", "Kendall", " Dcov", "Hoeffding", "MI", "IDA", "MIC", "RDC", "Lasso", "Elastic", "Z-score", "ProMISe")
	}
	
	if(is.matrix(TransResult)){
		TransResult=TransResult[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
		colnames(TransResult)=c( "Pearson", "Spearman", "Kendall", " Dcov", "Hoeffding", "MI", "IDA", "MIC", "RDC", "Lasso", "Elastic", "Z-score", "ProMISe")
	} else {
		tt=as.matrix(TransResult)
		TransResult=t(tt)
		TransResult=TransResult[,c(2,4,6,8,10,12,14,16,18,20,22,24,26)]
		#print(TransResult)
		names(TransResult)=c( "Pearson", "Spearman", "Kendall", " Dcov", "Hoeffding", "MI", "IDA", "MIC", "RDC", "Lasso", "Elastic", "Z-score", "ProMISe")
	}
	if(is.matrix(ExpResult)){
		numberExpResult=apply(ExpResult, 2, sum)
		ranking=apply(ExpResult, 1, function(i) rank(i))
		ranking=t(ranking)
		rankExpResult=apply(ranking, 2, sum)
	}
	if(is.matrix(TransResult)){
		numberTransResult=apply(TransResult, 2, sum)	
		rankingT=apply(TransResult, 1, function(i) rank(i))
		rankingT=t(rankingT)
		rankTransResult=apply(rankingT, 2, sum)
	}
	if(is.matrix(ExpResult)){
		Exp=list(ExpResult, numberExpResult, rankExpResult)
	} else Exp=ExpResult
	
	if(is.matrix(TransResult)){
		Trans=list(TransResult, numberTransResult, rankTransResult)
	} else Trans=TransResult
	
	result=list(Exp, Trans)
}

## Enrichment analysis including GO and KEGG enrichment analysis using two functions: GOBPenrichment and KEGGenrichment ##
# The input is a list of gene symbols. For example: Genes <- c("AREG", "FKBP5", "CXCL13", "KLF9", "ZC3H12A", "P4HA1", "TLE1", "CREB3L2", "TXNIP", "PBX1", "GJA1", "ITGB8", "CCL3", "CCND2", "KCNJ15", "CFLAR", "CXCL10", "CYSLTR1", "IGFBP7", "RHOB", "MAP3K5", "CAV2", "CAPN2", "AKAP13", "RND3", "IL6ST", "RGS1", "IRF4", "G3BP1", "SEL1L", "VEGFA", "SMAD1", "CCND1", "CLEC3B", "NEB", "AMD1", "PDCD4", "SCD", "TM2D3", "BACH2", "LDLR", "BMPR1B", "RFXAP", "ASPH", "PTK2B", "SLC1A5", "ENO2", "TRPM8", "SATB1", "MIER1", "SRSF1", "ATF3", "CCL5", "MCM6", "GCH1", "CAV1", "SLC20A1")
#@import GSEABase @import Category
## GO BP enrichment analysis, Cutoff is normally set to 0.05 ##
#' Functional enrichment analysis
#' 
#' GO BP enrichment analysis for a gene list
# @importMethodsFrom AnnotationDbi mget Lkeys get
# @importFrom org.Hs.eg.db org.Hs.egSYMBOL2EG org.Hs.egGO org.Hs.egSYMBOL
# @import GOstats 
# @import Category 
#' @param Genes a list of gene symbols
#' @param Cutoff the significant level, e.g. 0.05
#' @examples
#'  print("result = GOBPenrichment(genelist, 0.05)")
#' @return a list of GO terms for the genes
#' @export
#' @references
#' Ashburner, M., Ball, C.A., Blake, J.A., Botstein, D., Butler, H., Cherry, J.M., Davis, A.P., Dolinski, K., Dwight, S.S., Eppig, J.T., Harris, M.A., Hill, D.P., Issel-Tarver, L., Kasarskis, A., Lewis, S., Matese, J.C., Richardson, J.E., Ringwald, M., Rubin, G.M. and Sherlock, G. (2000) Gene Ontology: tool for the unification of biology. Nat. Genet., 25, 25-29.
GOBPenrichment <- function(Genes, Cutoff){
  
if(requireNamespace("AnnotationDbi", quitely=TRUE)) {
  EntrezIDs <- AnnotationDbi::mget(Genes, org.Hs.eg.db::org.Hs.egSYMBOL2EG, ifnotfound=NA)
  EntrezIDs <- as.character(EntrezIDs)
  GOAnnotation <- AnnotationDbi::get("org.Hs.egGO")
  Universe <- AnnotationDbi::Lkeys(GOAnnotation)

  Params <- new("GOHyperGParams",
                  geneIds=EntrezIDs,
                  universeGeneIds=Universe,
                  annotation="org.Hs.eg.db",
                  ontology="BP",
                  pvalueCutoff=Cutoff,
                  conditional=FALSE,
                  testDirection="over")
}
if(requireNamespace("GOstats", quitely=TRUE)) {
  Over <- GOstats::hyperGTest(Params)
}
if(requireNamespace("Category", quitely=TRUE)) {
  Glist <- Category::geneIdsByCategory(Over)
}
if(requireNamespace("AnnotationDbi", quitely=TRUE)){
Glist <- sapply(Glist, function(.ids) {
 	.sym <- AnnotationDbi::mget(.ids, envir=org.Hs.eg.db::org.Hs.egSYMBOL, ifnotfound=NA)
 	.sym[is.na(.sym)] <- .ids[is.na(.sym)]
 	paste(.sym, collapse=";")
 	})
BP <- summary(Over)
BP$Symbols <- Glist[as.character(BP$GOBPID)]
# Adjust p-value using Benjamini & Hochberg (BH) method
BP$adjustPvalue <-p.adjust(BP$Pvalue, "BH", length(BP$Pvalue))
# write.csv(BP,'BPResult.csv')
}

return(BP)
}

## KEGG enrichment analysis, Cutoff is normally set to 0.05 ##
## GO BP enrichment analysis, Cutoff is normally set to 0.05 ##
#' Functional enrichment analysis
#' KEGG enrichment analysis for a gene list
# @importMethodsFrom AnnotationDbi mget Lkeys get
# @import GOstats 
# @import Category 
# @importFrom org.Hs.eg.db org.Hs.egSYMBOL2EG org.Hs.egPATH org.Hs.egSYMBOL
#' @param Genes a list of gene symbols
#' @param Cutoff the significant level, e.g. 0.05
#' @examples
#'  print("result = KEGGenrichment(genelist, 0.05)") 
#' @return a list of pathways for the genes
#' @export
#' @references
#' Kanehisa, M. and Goto, S. (2000) KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res., 28, 27-30.

KEGGenrichment <- function(Genes, Cutoff){

  if(requireNamespace("AnnotationDbi", quitely=TRUE)){
EntrezIDs <- AnnotationDbi::mget(Genes, org.Hs.eg.db::org.Hs.egSYMBOL2EG, ifnotfound=NA)
EntrezIDs <- as.character(EntrezIDs)
KEGGAnnotation <- AnnotationDbi::get("org.Hs.egPATH")
Universe <- AnnotationDbi::Lkeys(KEGGAnnotation)
Params <- new("KEGGHyperGParams", 
                     geneIds=EntrezIDs, 
                     universeGeneIds=Universe, 
                     annotation="org.Hs.eg.db", 
                     categoryName="KEGG", 
                     pvalueCutoff=Cutoff,
                     testDirection="over")
}

if(requireNamespace("GOstats", quitely=TRUE)){
  Over <- GOstats::hyperGTest(Params)
 KEGG <- summary(Over)
}
if(requireNamespace("Category", quitely=TRUE)){
  Glist <- Category::geneIdsByCategory(Over)
}
if(requireNamespace("AnnotationDbi", quitely=TRUE)){
Glist <- sapply(Glist, function(.ids) {
 	.sym <- AnnotationDbi::mget(.ids, envir=org.Hs.eg.db::org.Hs.egSYMBOL, ifnotfound=NA)
 	.sym[is.na(.sym)] <- .ids[is.na(.sym)]
 	paste(.sym, collapse=";")
 	})
KEGG$Symbols <- Glist[as.character(KEGG$KEGGID)]
# Adjust p-value using Benjamini & Hochberg (BH) method
KEGG$adjustPvalue <-p.adjust(KEGG$Pvalue, "BH", length(KEGG$Pvalue))
# write.csv(KEGG,'KEGGResult.csv')
}

return(KEGG)
}

#' Function for validate the results from all 12 methods.
#' @param allmethods A list of results (matrix with columns are miRNA and rows are mRNAs). 
#' @param topk Top k targets of each miRNA that will be extracted for validation
#' @param Expgroundtruth The ground truth in .csv file for validation
#' @param LFC log fold-change for validating the results using transfection experiments
#' @param downreg If set to TRUE the negative effects will have higher ranks than the positives.
#' @return The validation results for all 12 methods 
#' @export
experiment=function(allmethods, topk, Expgroundtruth, LFC, downreg){
  
  psv=ValidateAll(allmethods[[1]], topk, Expgroundtruth, LFC, downreg)
  spv=ValidateAll(allmethods[[2]], topk, Expgroundtruth, LFC, downreg)
  kendallv=ValidateAll(allmethods[[3]], topk, Expgroundtruth, LFC, downreg)
  dcovv=ValidateAll(allmethods[[4]], topk, Expgroundtruth, LFC, downreg)
  hoeffdingv=ValidateAll(allmethods[[5]], topk, Expgroundtruth, LFC, downreg)
  miv=ValidateAll(allmethods[[6]], topk, Expgroundtruth, LFC, downreg)
  idav=ValidateAll(allmethods[[7]], topk, Expgroundtruth, LFC, downreg)
 # micv=ValidateAll(allmethods[[8]], topk, Expgroundtruth, LFC, downreg, TargetBinding)
  rdcv=ValidateAll(allmethods[[8]], topk, Expgroundtruth, LFC, downreg)
  lassov=ValidateAll(allmethods[[9]], topk, Expgroundtruth, LFC, downreg)
  elasticv=ValidateAll(allmethods[[10]], topk, Expgroundtruth, LFC, downreg)
  zsv=ValidateAll(allmethods[[11]], topk, Expgroundtruth, LFC, downreg)
  promisev=ValidateAll(allmethods[[12]], topk, Expgroundtruth, LFC, downreg)
  
  #############genenames
  miRs=colnames(allmethods[[1]])
  #########decorate and return
  result1= cbind(psv[[2]], spv[[2]][,2:3], kendallv[[2]][,2:3], dcovv[[2]][,2:3], hoeffdingv[[2]][,2:3], miv[[2]][,2:3], idav[[2]][,2:3], rdcv[[2]][,2:3], lassov[[2]][,2:3], elasticv[[2]][,2:3],zsv[[2]][,2:3], promisev[[2]][,2:3])
  rownames(result1)=miRs
  result2= cbind(psv[[4]], spv[[4]][,2:3], kendallv[[4]][,2:3], dcovv[[4]][,2:3], hoeffdingv[[4]][,2:3], miv[[4]][,2:3], idav[[4]][,2:3], rdcv[[4]][,2:3], lassov[[4]][,2:3], elasticv[[4]][,2:3],zsv[[4]][,2:3], promisev[[4]][,2:3])
  rownames(result2)=miRs
  result=list(result1, result2)
  return(result)
}

## Compare the validation results of 13 built-in methods ##
#' Filter and compare the validation results from 12 methods
#' Keep the miRNAs that have at least noVal confirmed targets and compare the validation results from all methods.
#' @param allresults the results from all methods generated from experiment function. This is a list.
#' @param noVal Number of confirmed targets in each method (threshold) to filter. Records (miRNA) with less than this will be removed
#' @return the validation results of all methods
#' @examples
#' print("result=filterAndCompare(allresults, 2)")
#' @export 
filterAndCompare=function(allresults, noVal){
  #allresults: the results from all methods generated from experiment function. This is a list
  #noVal: number of confirmed targets in each method (threshold) to filter. Records (miRNA) with less than this will be removed
  ExpResult=allresults[[1]]
  TransResult=allresults[[2]]
  temp1=apply(ExpResult, 1, function(x) all(x[c(2,4,6,8,10,12,14,16,18,20,22,24)]>(noVal-1)))
  ExpResult=ExpResult[temp1,]
  temp2=apply(TransResult, 1, function(x) all(x[c(2,4,6,8,10,12,14,16,18,20,22,24)]>(noVal-1)))
  TransResult=TransResult[temp2,]
  ###If else to incase only one record. In that case the result is not a matrix
  if(is.matrix(ExpResult)){
    ExpResult=ExpResult[,c(2,4,6,8,10,12,14,16,18,20,22,24)]
    colnames(ExpResult)=c( "Pearson", "Spearman", "Kendall", " Dcov", "Hoeffding", "MI", "IDA", "RDC", "Lasso", "Elastic", "Z-score", "ProMISe")
  } else {
    tt=as.matrix(ExpResult)
    ExpResult=t(tt)
    ExpResult=ExpResult[,c(2,4,6,8,10,12,14,16,18,20,22,24)]
    names(ExpResult)=c( "Pearson", "Spearman", "Kendall", " Dcov", "Hoeffding", "MI", "IDA", "RDC", "Lasso", "Elastic", "Z-score", "ProMISe")
  }
  
  if(is.matrix(TransResult)){
    TransResult=TransResult[,c(2,4,6,8,10,12,14,16,18,20,22,24)]
    colnames(TransResult)=c( "Pearson", "Spearman", "Kendall", " Dcov", "Hoeffding", "MI", "IDA", "RDC", "Lasso", "Elastic", "Z-score", "ProMISe")
  } else {
    tt=as.matrix(TransResult)
    TransResult=t(tt)
    TransResult=TransResult[,c(2,4,6,8,10,12,14,16,18,20,22,24)]
    #print(TransResult)
    names(TransResult)=c( "Pearson", "Spearman", "Kendall", " Dcov", "Hoeffding", "MI", "IDA", "RDC", "Lasso", "Elastic", "Z-score", "ProMISe")
  }
  if(is.matrix(ExpResult)){
    numberExpResult=apply(ExpResult, 2, sum)
    ranking=apply(ExpResult, 1, function(i) rank(i))
    ranking=t(ranking)
    rankExpResult=apply(ranking, 2, sum)
  }
  if(is.matrix(TransResult)){
    numberTransResult=apply(TransResult, 2, sum)	
    rankingT=apply(TransResult, 1, function(i) rank(i))
    rankingT=t(rankingT)
    rankTransResult=apply(rankingT, 2, sum)
  }
  if(is.matrix(ExpResult)){
    Exp=list(ExpResult, numberExpResult, rankExpResult)
  } else Exp=ExpResult
  
  if(is.matrix(TransResult)){
    Trans=list(TransResult, numberTransResult, rankTransResult)
  } else Trans=TransResult
  
  result=list(Exp, Trans)
}

#' Convert miRNA symbols from a miRBase version to another 
#' 
#' This function convert the miRNAs in the input file from the "source" miRBase version to the "Target" version. 
#' If users do not know the miRBase version of the input file, please set the source version to 0. The function will match the 
#' miRNAs in the input file to all miRBase versions to find the most likely miRBase version. Currently, we have versions 16-21.
#' @param miRNAListFile the input file containing a list of miRNA symbols in csv format
#' @param sourceV the miRBase version of the input miRNAs, e.g. 16. If users do not know the version, use 0.
#' @param targetV the miRBase version that we want to convert into, e.g. 21.
#' @return A csv file in the working directory containing the converted miRNA symbols.
#' @examples 
#' miRs=system.file("extdata", "ToymiRs.csv", package="miRLAB")
#' convert(miRs, 17, 21) 
#' @export 
convert = function (miRNAListFile,sourceV,targetV) {
  load(system.file("extdata", "database.RData", package="miRLAB"))
  miRNAList = as.matrix( read.csv( miRNAListFile ,header = FALSE) )
  sourceName = miRNAList
  sourceVersion = c()
  targetName = c()
  targetVersion = c()
  
  if (sourceV != 0) # have the source version
  {
    location = match( miRNAList, all[,sourceV-14] )
    isNA = is.na(location)
    targetVersion[which(!isNA)] = targetV
    targetVersion[which(isNA)] = NA
    targetName = all[location, targetV-14]
    sourceVersion = rep(sourceV,length(miRNAList))
  }else 
  {
    allVersionList = c(all[,2],all[,3],all[,4],all[,5],all[,6],all[,7])
    location = match(miRNAList, allVersionList)
    sourceVersion = 16 + (location %/% 2602)
    isNA = is.na(location)
    targetVersion[which(!isNA)] = targetV
    targetVersion[which(isNA)] = NA
    location = location %% 2602
    location[which(location == 0)] = 2602
    targetName = all[location, targetV-14]
  }
  
  res = cbind(sourceName, sourceVersion, targetName, targetVersion)
  colnames(res) = c("sourceName","sourceVersion","targetName","targetVersion")
  write.table(res, file="resOfConvert.csv", sep=",",row.names = FALSE,col.names = TRUE)
}

# A function for preparing data to classify cancer subtypes
#' @importFrom utils write.table
prepareDataForPam50 = function(d, dataDir = "bioclassifier_data", inputFileName = "inputFile.txt") {
  
  d <- t(d)
  noOfRow <- nrow(d)
  noOfCol <- ncol(d)
  temp <- matrix(0, nrow = (noOfRow + 1), ncol = noOfCol)
  row.names(temp) <- c("T", row.names(d))
  colnames(temp) <- colnames(d)
  temp[2:(noOfRow + 1), 1:noOfCol] <- d[, ]
  dir.create(file.path(".", dataDir), showWarnings = FALSE)
  write.table(temp, paste("./", dataDir, "/", inputFileName, sep = ""), col.names=NA, sep = "\t")
}

#' Identify miRNA targets by ICP and PAM50
#' 
#' This function identifies miRNA targets by ICP and PAM50.
#' @importFrom utils read.table
#' @param d A matrix of expression of miRNAs and mRNAs with columns being miRNA or mRNA names and rows being samples
#' @param nmiR Number of miRNAs
#' @param nmR Number of mRNAs
#' @param fiftymRNAsData A matrix of expression of 50 mRNAs in PAM50 with columns being mRNA names and rows being samples
#' @return The matrix of causal effects of miRNAs and mRNAs with columns being miRNAs and rows being mRNAs
#' @references 
#' 1. Parker, J. S., et al. (2009). "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes." Journal of Clinical Oncology 27(8): 1160-1167.
#' @export 
ICPPam50 = function(d, nmiR, nmR, fiftymRNAsData) {
  
  # Set parameters
  ## Folder including input and output files used to categorise samples into cancer subtypes using Pam50
  dataDir = "bioclassifier_data"
  ## Input file used to categorise samples into cancer subtypes using Pam50
  inputFileName = "inputFile.txt"
  
  # Classify samples based on cancer subtypes using Pam50
  # Reference: Parker, J. S., et al. (2009). "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes." Journal of Clinical Oncology 27(8): 1160-1167.
  ## Prepare data
  prepareDataForPam50(fiftymRNAsData, dataDir, inputFileName)
  ## Run the assignment algorithm
  source(system.file(package = 'miRLAB', 'bioclassifier_R', 'subtypePrediction_functions.R'))
  source(system.file(package = 'miRLAB', 'bioclassifier_R', 'subtypePrediction_distributed.R'))
  
  ## Read the result
  result <- read.table(paste("./", dataDir, "/", "outputFile_pam50scores.txt", sep = ""))
  result <- result[-1, -1]
  
  # Define experimental settings
  ExpInd <- c(rep(0, nrow(d)))
  ExpInd[which(result[, 6] == "Basal")] <- 1
  ExpInd[which(result[, 6] == "Her2")] <- 2
  ExpInd[which(result[, 6] == "LumA")] <- 3
  ExpInd[which(result[, 6] == "LumB")] <- 4
  ExpInd[which(result[, 6] == "Normal")] <- 5
  
  # Estimate causal effects
  standardizedData <- scale(d)
  temp = matrix(nrow = nmR, ncol = nmiR)
  for (i in 1:nmR) {
    Y <- standardizedData[, nmiR+i]
    X <- standardizedData[, c(1:nmiR)]
    icp <- InvariantCausalPrediction::hiddenICP(X, Y, ExpInd, alpha = 0.1, mode = "asymptotic", intercept=FALSE)
    for (k in 1:nmiR) {
      temp[i,k] <- icp$betahat[k]
    }
  }
  colnames(temp) <- colnames(d)[1:nmiR]
  row.names(temp) <- colnames(d)[(nmiR + 1):(nmiR + nmR)]
  
  return(temp)
}

#' Identify the top miRNA targets by ICP and PAM50
#' 
#' This function identifies the top miRNA targets by ICP and PAM50.
#' @param d A matrix of expression of miRNAs and mRNAs with columns being miRNA or mRNA names and rows being samples
#' @param nmiR Number of miRNAs
#' @param nmR Number of mRNAs
#' @param fiftymRNAsData A matrix of expression of 50 mRNAs in PAM50 with columns being mRNA names and rows being samples
#' @param top 1 if getting the top of all miRNAs and 2 if getting the top of each miRNA
#' @param topk Number of the top to get
#' @return The top k miRNA targets
#' @references 
#' 1. Parker, J. S., et al. (2009). "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes." Journal of Clinical Oncology 27(8): 1160-1167.
#' @export 
identifymiRTargetsByICPPam50 = function(d, nmiR, nmR, fiftymRNAsData, top = 1, topk = 500) {
  
  temp <- ICPPam50(d, nmiR, nmR, fiftymRNAsData)
  
  # Get the result
  miRTargets = NULL
  if(top == 1) {# Top for all miRNAs
    miRTargets = Extopk(temp, topk)
  } else { # Top for each miRNA
    for(i in 1:nmiR) {
      miRiTopk = bRank(temp, i, topk, FALSE)
      miRTargets <- rbind(miRTargets, miRiTopk)
    }
  }
  
  return(miRTargets)
}

#' Identify the top miRNA targets by an ensemble method with ICP-PAM50, Pearson and Lasso
#' 
#' This function identifies the top miRNA targets by an ensemble method with ICP-PAM50, Pearson and Lasso.
#' @param d A matrix of expression of miRNAs and mRNAs with columns being miRNA or mRNA names and rows being samples
#' @param nmiR Number of miRNAs
#' @param nmR Number of mRNAs
#' @param fiftymRNAsData A matrix of expression of 50 mRNAs in PAM50 with columns being mRNA names and rows being samples
#' @param top 1 if getting the top of all miRNAs and 2 if getting the top of each miRNA
#' @param topk Number of the top to get
#' @return The top k miRNA targets
#' @references 
#' 1. Parker, J. S., et al. (2009). "Supervised Risk Predictor of Breast Cancer Based on Intrinsic Subtypes." Journal of Clinical Oncology 27(8): 1160-1167.
#' @export 
identifymiRTargetsByEnsemble = function(d, nmiR, nmR, fiftymRNAsData, top = 1, topk = 500) {
  
  # Get the results of ICPPam50, Pearson and Lasso
  hid <- ICPPam50(d, nmiR, nmR, fiftymRNAsData)
  standardizedData <- scale(d)
  dir.create(file.path(".", "temp_data"), showWarnings = FALSE)
  write.table(standardizedData, file = "./temp_data/input.csv", sep = ",", row.names = FALSE)
  pea = Pearson("./temp_data/input.csv", cause = 1:nmiR, effect = (nmiR+1):(nmiR+nmR))
  las = Lasso("./temp_data/input.csv", cause = 1:nmiR, effect = (nmiR+1):(nmiR+nmR))
  
  row.names(pea) <- row.names(hid)
  row.names(las) <- row.names(hid)
  
  # Call the ensemble method
  borda = Borda(list(hid, pea, las))
  
  # Get the result
  miRTargets = NULL
  if(top == 1) {# Top for all miRNAs
    miRTargets = Extopk(borda, topk)
  } else { # Top for each miRNA
    for(i in 1:nmiR) {
      miRiTopk = bRank(borda, i, topk, FALSE)
      miRTargets <- rbind(miRTargets, miRiTopk)
    }
  }
  
  return(miRTargets)
}

				