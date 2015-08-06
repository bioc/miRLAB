# This file contains a collection of methods for inferring miRNA-mRNA regulatory relationships
# Each method will take the dataset (.csv) as an input with the first n colomns are miRNAs 
# and the last m columns are mRNAs. Rows are samples. We need to specify the number of miRNAs in the dataset
# and the index of the miRNA, i.e. column that contains the miNRA. Topk is the number of top k interactions we
# would like to extract.

###########################################################TCGAAssembler- Module_A.R#####################################
#############################################################################

# TCGA-Assembler : An open-source R program for downloading, processing and analyzing public TCGA data.
# Copyright (C) <2014>  <Yitan Zhu>
# This file is part of TCGA-Assembler.

#  TCGA-Assembler is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.

#  TCGA-Assembler is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#  You should have received a copy of the GNU General Public License
#  along with TCGA-Assembler.  If not, see <http://www.gnu.org/licenses/>.
############################################################################  


##############################################################################

# TCGA-Assembler Version 1.0.3 Module A

##############################################################################

#' @import RCurl
#' @import httr
#' @import stringr





DownloadRNASeqData <-function (traverseResultFile, saveFolderName, cancerType, assayPlatform, dataType = "", tissueType = NULL, inputPatientIDs = NULL, outputFileName = "")
{
  
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download Gene Expression data of ", cancerType, " patients.", sep = ""));
  
  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }
  
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }  
  
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  file_url=c() 
  dir_url=c()
  upper_file_url=c()
  load(traverseResultFile);  
  dir.create(path = saveFolderName, recursive = TRUE);
  if (assayPlatform == "RNASeqV1")
  {
    platform = c("illuminaga_rnaseq", "illuminahiseq_rnaseq"); 
    Institution = c("unc.edu", "bcgsc.ca");    
  }
  if (assayPlatform == "RNASeqV2")
  {
    platform = c("illuminaga_rnaseqv2", "illuminahiseq_rnaseqv2"); 
    Institution = c("unc.edu");    
  }
  if (assayPlatform == "Microarray")
  {
    platform = c("agilentg4502a_07_3", "ht_hg-u133a", "agilentg4502a_07_1", "agilentg4502a_07_2", "hg-u133_plus_2");
    Institution = c("unc.edu", "broad.mit.edu", "genome.wustl.edu");
    dataType = "";
  }
  
  # For returning downloaded data
  downloadedData = vector("list", 0);     
  dataIndex = 0;     
  
  # download RNASeqV2 data
  if (assayPlatform == "RNASeqV2")
  {
    for (IDin in 1:length(Institution))
    {
      for (IDpl in 1:length(platform))
      {
        SpecificName = paste(cancerType, "__", Institution[IDin], "__", platform[IDpl], sep = "");
        SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/", Institution[IDin], "/", platform[IDpl], "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
        RNALevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        sdrf = toupper(downloadResult$data);
        
        level_3_filename_column = max(grep(pattern = "Derived Data File", x = sdrf[1, ], ignore.case = TRUE));
        DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
        ExtractNameColID = grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE);
        RefGenomeColID = grep(pattern = "Comment \\[Genome reference]", x = sdrf[1, ], ignore.case = TRUE);
        if (length(ExtractNameColID) == 0)
        {
          ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));          
        }
        colnames(sdrf) = sdrf[1, ];
        sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
        Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
        if (length(Level3_ID) == 0)
        {
          next;
        }
        sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column, RefGenomeColID), drop = FALSE];
        sdrf = sdrf[!duplicated(sdrf[, 2]), , drop = FALSE];
        
        # Only keep the file information for the data types that should be downloaded.
        keepID = c();
        for (keep_i in 1:length(dataType))
        {
          keepID = c(keepID, grep(pattern = dataType[keep_i], x = sdrf[, 2], ignore.case = TRUE))
        }
        sdrf = sdrf[sort(unique(keepID), decreasing = FALSE), , drop = FALSE];
        
        # If specific patient TCGA barcodes are inputted, only download the specified samples.
        if (!is.null(inputPatientIDs))
        {
          indInputPatientID = c();
          for (i in 1:length(inputPatientIDs))
          {
            indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
          }
          if (length(indInputPatientID) == 0)
          {
            next;
          }else{
            sdrf = sdrf[indInputPatientID, , drop = FALSE];
          }
        }
        
        # Download data of specified tissue
        if (!is.null(tissueType))
        {
          SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                             Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
          sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
        }
        
        if (dim(sdrf)[1] == 0)
        {
          next;
        }
        
        # Start to download data files
        exp_gene = NULL;
        column_gene = NULL;
        gene_left_column = NULL;
        gene_data = NULL;
        exp_gene_normalized = NULL;
        column_gene_normalized = NULL;
        gene_normalized_left_column = NULL;
        gene_normalized_data = NULL;    
        
        exp_isoform = NULL;
        column_isoform = NULL;
        isoform_left_column = NULL;
        isoform_data = NULL;
        exp_isoform_normalized = NULL;
        column_isoform_normalized = NULL;
        isoform_normalized_left_column = NULL;
        isoform_normalized_data = NULL;    
        
        exp_exon = NULL;
        column_exon = NULL;
        exon_left_column = NULL;
        exon_data=NULL;
        exp_junction = NULL;
        column_junction = NULL;
        junction_left_column = NULL;
        junction_data=NULL;
        
        for (i in 1:dim(sdrf)[1])
        {
          time1 = proc.time();
          
          sample_TCGA_id = sdrf[i, 1];
          ind = RNALevel3ID[grepEnd(pattern = sdrf[i, 2], x = upper_file_url[RNALevel3ID], ignore.case = FALSE)];
          if (length(ind) == 0)
          {
            next;
          }
          if (length(ind) > 1)
          {
            URL = GetNewestURL(AllURL = file_url[ind]);
          }else{
            URL = file_url[ind];
          }  
          downloadResult = urlReadTable(url = URL);
          if (downloadResult$errorFlag != 0)
          {
            next;
          }
          s = downloadResult$data;
          
          # read gene expression data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("rsem.genes.results"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(gene_left_column))
              {
                gene_left_column = s[, c(1, 4), drop = FALSE];
                gene_data = s[, c(2, 3), drop = FALSE];
                exp_gene = c(sample_TCGA_id, sample_TCGA_id);
                column_gene = c("raw_count", "scaled_estimate");
              }else{
                if (sum(gene_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                gene_data = cbind(gene_data, s[, c(2, 3), drop = FALSE]);
                exp_gene = c(exp_gene, sample_TCGA_id, sample_TCGA_id);
                column_gene = c(column_gene, "raw_count", "scaled_estimate");
              }
            }
          }
          
          # read gene normalized data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("rsem.genes.normalized_results"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(gene_normalized_left_column))
              {
                gene_normalized_left_column = s[, 1, drop = FALSE];
                gene_normalized_data = s[, 2, drop = FALSE];
                exp_gene_normalized = sample_TCGA_id;
                column_gene_normalized = "normalized_count";
              }else{
                if (sum(gene_normalized_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                gene_normalized_data = cbind(gene_normalized_data, s[, 2, drop = FALSE]);
                exp_gene_normalized = c(exp_gene_normalized, sample_TCGA_id);
                column_gene_normalized = c(column_gene_normalized, "normalized_count");
              }
            }
          }          
          
          # read isoform data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("rsem.isoforms.results"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(isoform_left_column))
              {
                isoform_left_column = s[, 1, drop = FALSE];
                isoform_data = s[, c(2, 3), drop = FALSE];
                exp_isoform = c(sample_TCGA_id, sample_TCGA_id);
                column_isoform = c("raw_count", "scaled_estimate");
              }else{
                if (sum(isoform_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                isoform_data = cbind(isoform_data, s[, c(2, 3), drop = FALSE]);
                exp_isoform = c(exp_isoform, sample_TCGA_id, sample_TCGA_id);
                column_isoform = c(column_isoform, "raw_count", "scaled_estimate");
              }
            }
          }
          
          # read isoform normalized data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("rsem.isoforms.normalized_results"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(isoform_normalized_left_column))
              {
                isoform_normalized_left_column = s[, 1, drop = FALSE];
                isoform_normalized_data = s[, 2, drop = FALSE];
                exp_isoform_normalized = sample_TCGA_id;
                column_isoform_normalized = "normalized_count";
              }else{
                if (sum(isoform_normalized_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                isoform_normalized_data = cbind(isoform_normalized_data, s[, 2, drop = FALSE]);
                exp_isoform_normalized = c(exp_isoform_normalized, sample_TCGA_id);
                column_isoform_normalized = c(column_isoform_normalized, "normalized_count");
              }
            }
          }          
          
          # read exon data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("exon_quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(exon_left_column))
              {
                exon_left_column = s[, 1, drop = FALSE];
                exon_data = s[, 2:4, drop = FALSE];
                exp_exon = c(sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_exon = c("raw_counts", "median_length_normalized", "RPKM");
              }else{
                if (sum(exon_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                exon_data = cbind(exon_data, s[, 2:4, drop = FALSE]);
                exp_exon = c(exp_exon, sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_exon = c(column_exon, "raw_counts", "median_length_normalized", "RPKM");
              }
            }
          }
          
          # read junction data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("junction_quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(junction_left_column))
              {
                junction_left_column = s[, 1, drop = FALSE];
                junction_data = s[, 2, drop = FALSE];
                exp_junction = sample_TCGA_id;
                column_junction = "raw_counts";
              }else{
                if (sum(junction_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                junction_data = cbind(junction_data, s[, 2, drop = FALSE]);
                exp_junction = c(exp_junction, sample_TCGA_id);
                column_junction = c(column_junction, "raw_counts");
              }
            }
          }          
          
          time = proc.time() - time1;
          writeLines(paste("Downloaded - ", SpecificName, " - file ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        
        writeLines("Save data to local disk.");
        ID = str_locate_all(traverseResultFile, "_")[[1]];
        ID = ID[dim(ID)[1], 2];
        if (length(grep(pattern = "rsem.genes.results", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__rsem.genes.results__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", "Hybridization REF", exp_gene), c("gene_id", "transcript_id", column_gene), cbind(gene_left_column, gene_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE); 
          
          # For returning downloaded data
          dataIndex = dataIndex + 1
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", "Hybridization REF", exp_gene), c("gene_id", "transcript_id", column_gene), cbind(gene_left_column, gene_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "rsem.genes.results", sep = "__");
        }
        if (length(grep(pattern = "rsem.genes.normalized_results", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__rsem.genes.normalized_results__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_gene_normalized), c("gene_id", column_gene_normalized), cbind(gene_normalized_left_column, gene_normalized_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE); 
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_gene_normalized), c("gene_id", column_gene_normalized), cbind(gene_normalized_left_column, gene_normalized_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "rsem.genes.normalized_results", sep = "__");
        }
        if (length(grep(pattern = "rsem.isoforms.results", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__rsem.isoforms.results__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_isoform), c("isoform_id", column_isoform), cbind(isoform_left_column, isoform_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE); 
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_isoform), c("isoform_id", column_isoform), cbind(isoform_left_column, isoform_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "rsem.isoforms.results", sep = "__");
        }
        if (length(grep(pattern = "rsem.isoforms.normalized_results", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__rsem.isoforms.normalized_results__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_isoform_normalized), c("isoform_id", column_isoform_normalized), cbind(isoform_normalized_left_column, isoform_normalized_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_isoform_normalized), c("isoform_id", column_isoform_normalized), cbind(isoform_normalized_left_column, isoform_normalized_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "rsem.isoforms.normalized_results", sep = "__");
        }
        if (length(grep(pattern = "exon_quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__exon_quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_exon), c("exon", column_exon), cbind(exon_left_column, exon_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_exon), c("exon", column_exon), cbind(exon_left_column, exon_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "exon_quantification", sep = "__");
        }
        if (length(grep(pattern = "junction_quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__junction_quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_junction), c("junction", column_junction), cbind(junction_left_column, junction_data)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE); 
          
          # For returning downloaded data
          dataIndex = dataIndex + 1          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_junction), c("junction", column_junction), cbind(junction_left_column, junction_data));
          names(downloadedData)[dataIndex] = paste(SpecificName, "junction_quantification", sep = "__");
        }
      }
    }
  }
  
  # download RNASeqV1 data
  if (assayPlatform == "RNASeqV1")
  {
    for (IDin in 1:length(Institution))
    {
      for (IDpl in 1:length(platform))
      {
        SpecificName = paste(cancerType, "__", Institution[IDin], "__", platform[IDpl], sep = "");
        SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/", Institution[IDin], "/", platform[IDpl], "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
        RNALevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        sdrf = toupper(downloadResult$data);
        
        level_3_filename_column = max(grep(pattern = "Derived Data File", x = sdrf[1, ], ignore.case = TRUE));
        DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
        ExtractNameColID = grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE);
        RefGenomeColID = grep(pattern = "Comment \\[Genome reference]", x = sdrf[1, ], ignore.case = TRUE);
        if (length(ExtractNameColID) == 0)
        {
          ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));          
        }
        colnames(sdrf) = sdrf[1, ];
        sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
        Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
        if (length(Level3_ID) == 0)
        {
          next;
        }
        sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column, RefGenomeColID), drop = FALSE];
        sdrf = sdrf[!duplicated(sdrf[, 2]), , drop = FALSE];
        
        # Only keep the file information for the data types that should be downloaded.
        keepID = c();
        for (keep_i in 1:length(dataType))
        {
          keepID = c(keepID, grep(pattern = dataType[keep_i], x = sdrf[, 2], ignore.case = TRUE))
        }
        sdrf = sdrf[sort(unique(keepID), decreasing = FALSE), , drop = FALSE];
        
        # If specific patient TCGA barcodes are inputted, only download the specified samples.
        if (!is.null(inputPatientIDs))
        {
          indInputPatientID = c();
          for (i in 1:length(inputPatientIDs))
          {
            indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
          }
          if (length(indInputPatientID) == 0)
          {
            next;
          }else{
            sdrf = sdrf[indInputPatientID, , drop = FALSE];
          }
        }
        
        # Download data of specified tissue
        if (!is.null(tissueType))
        {
          SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                             Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
          sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
        }
        
        if (dim(sdrf)[1] == 0)
        {
          next;
        }
        
        exp_names_gene = NULL;
        column_gene = NULL;
        gene_left_column = NULL;
        gene_RPKM = NULL;
        
        exp_names_exon = NULL;
        column_exon = NULL;
        exon_left_column = NULL;
        exon_RPKM = NULL;
        
        junction_count = NULL;
        exp_names_junction = NULL;    
        column_junction = NULL;     
        junction_left_column = NULL;        
        
        for (i in 1:dim(sdrf)[1])
        {
          time1 = proc.time();
          
          sample_TCGA_id = sdrf[i, 1];
          ind = RNALevel3ID[grepEnd(pattern = sdrf[i, 2], x = upper_file_url[RNALevel3ID], ignore.case = FALSE)];
          if (length(ind) == 0)
          {
            next;
          }
          if (length(ind) > 1)
          {
            URL = GetNewestURL(AllURL = file_url[ind]);
          }else{
            URL = file_url[ind];
          }  
          downloadResult = urlReadTable(url = URL);
          if (downloadResult$errorFlag != 0)
          {
            next;
          }
          s = downloadResult$data;
          
          # read gene expression data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("gene.quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(gene_left_column))
              {
                gene_left_column = s[, 1, drop = FALSE];
                gene_RPKM = s[, 2:4, drop = FALSE];
                exp_names_gene = c(sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_gene = c("raw_counts", "median_length_normalized", "RPKM");
              }else{
                if (sum(gene_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                gene_RPKM = cbind(gene_RPKM, s[, 2:4, drop = FALSE]);
                exp_names_gene = c(exp_names_gene, sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_gene = c(column_gene, "raw_counts", "median_length_normalized", "RPKM");
              }
            }
          }
          
          # read exon expression data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("exon.quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(exon_left_column))
              {
                exon_left_column = s[, 1, drop = FALSE];
                exon_RPKM = s[, 2:4, drop = FALSE];
                exp_names_exon = c(sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_exon = c("raw_counts", "median_length_normalized", "RPKM");
              }else{
                if (sum(exon_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                exon_RPKM = cbind(exon_RPKM, s[, 2:4, drop = FALSE]);
                exp_names_exon = c(exp_names_exon, sample_TCGA_id, sample_TCGA_id, sample_TCGA_id);
                column_exon = c(column_exon, "raw_counts", "median_length_normalized", "RPKM");
              }
            }
          }
          
          # read junction expression data
          for (jj in 1:1)
          {
            if (length(grep(pattern = toupper("spljxn.quantification"), x = sdrf[i, 2], ignore.case = FALSE)) > 0)
            {
              s = s[2:dim(s)[1], , drop = FALSE];
              I_order_probes = order(s[, 1], decreasing = FALSE);
              s = s[I_order_probes, , drop = FALSE];
              if (is.null(junction_left_column))
              {
                junction_left_column = s[, 1, drop = FALSE];
                junction_count = s[, 2, drop = FALSE];
                exp_names_junction = sample_TCGA_id;
                column_junction = "raw_counts";
              }else{
                if (sum(junction_left_column[, 1] != s[, 1]) > 0)
                {
                  next;
                }
                junction_count = cbind(junction_count, s[, 2, drop = FALSE]);
                exp_names_junction = c(exp_names_junction, sample_TCGA_id);
                column_junction = c(column_junction, "raw_counts");
              }
            }
          }
          
          time = proc.time() - time1;
          writeLines(paste("Downloaded - ", SpecificName, " - file ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        
        writeLines("Save data to local disk.");
        ID = str_locate_all(traverseResultFile, "_")[[1]];
        ID = ID[dim(ID)[1], 2];
        if (length(grep(pattern = "gene.quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__gene.quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names_gene), c("gene", column_gene), cbind(gene_left_column, gene_RPKM)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex + 1;          
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names_gene), c("gene", column_gene), cbind(gene_left_column, gene_RPKM));
          names(downloadedData)[dataIndex] = paste(SpecificName, "gene.quantification", sep = "__");
        }
        if (length(grep(pattern = "exon.quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__exon.quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names_exon), c("exon", column_exon), cbind(exon_left_column, exon_RPKM)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex + 1;
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names_exon), c("exon", column_exon), cbind(exon_left_column, exon_RPKM));
          names(downloadedData)[dataIndex] = paste(SpecificName, "exon.quantification", sep = "__");          
        }
        if (length(grep(pattern = "spljxn.quantification", x = dataType, ignore.case = TRUE)) > 0)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__spljxn.quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names_junction), c("junction", column_junction), cbind(junction_left_column, junction_count)), 
                      file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex + 1;
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names_junction), c("junction", column_junction), cbind(junction_left_column, junction_count));
          names(downloadedData)[dataIndex] = paste(SpecificName, "spljxn.quantification", sep = "__");             
        }
      }
    }
  }
  
  # download Microarray data
  if (assayPlatform == "Microarray")
  {
    for (IDin in 1:length(Institution))
    {
      for (IDpl in 1:length(platform))
      {
        
        SpecificName = paste(cancerType, "__", Institution[IDin], "__", platform[IDpl], sep = "");
        SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/", Institution[IDin], "/", platform[IDpl], "/", sep = "")), x = upper_file_url, ignore.case = FALSE);
        RNALevel3ID = SpecificID[grep(pattern = toupper("Level_3"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        sdrf = toupper(downloadResult$data);
        
        # Derived Array Data Matrix File
        level_3_filename_column = max(grep(pattern = "Derived Array Data Matrix File", x = sdrf[1, ], ignore.case = TRUE));
        DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
        ExtractNameColID = grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE);
        if (length(ExtractNameColID) == 0)
        {
          ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));          
        }
        colnames(sdrf) = sdrf[1, ];
        sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
        Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
        if (length(Level3_ID) == 0)
        {
          next;
        }
        sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column), drop = FALSE];
        sdrf = sdrf[!duplicated(sdrf[, 2]), , drop = FALSE];
        
        # If specific patient TCGA barcodes are inputted, only download the specified samples.
        if (!is.null(inputPatientIDs))
        {
          indInputPatientID = c();
          for (i in 1:length(inputPatientIDs))
          {
            indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
          }
          if (length(indInputPatientID) == 0)
          {
            next;
          }else{
            sdrf = sdrf[indInputPatientID, , drop = FALSE];
          }
        }
        
        # Download data of specified tissue
        if (!is.null(tissueType))
        {
          SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                             Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
          sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
        }
        
        if (dim(sdrf)[1] == 0)
        {
          next;
        }
        
        exp_names_gene = NULL;
        column_gene = NULL;
        gene_left_column = NULL;
        gene_RPKM = NULL;
        for (i in 1:dim(sdrf)[1])
        {
          time1 = proc.time();
          
          sample_TCGA_id = sdrf[i, 1];
          ind = RNALevel3ID[grepEnd(pattern = sdrf[i, 2], x = upper_file_url[RNALevel3ID], ignore.case = FALSE)];
          if (length(ind) == 0)
          {
            next;
          }
          if (length(ind) > 1)
          {
            URL = GetNewestURL(AllURL = file_url[ind]);
          }else{
            URL = file_url[ind];
          }  
          downloadResult = urlReadTable(url = URL);
          if (downloadResult$errorFlag != 0)
          {
            next;
          }
          s = downloadResult$data;
          
          column_gene = c(column_gene, s[2, 2]);
          s = s[3:dim(s)[1], , drop = FALSE];
          I_order_probes = order(s[, 1], decreasing = FALSE);
          s = s[I_order_probes, , drop = FALSE];
          if (is.null(gene_left_column))
          {
            gene_left_column = s[, 1, drop = FALSE];
            gene_RPKM = s[, 2, drop = FALSE];
            exp_names_gene = c(sample_TCGA_id);
          }else{
            if (sum(gene_left_column[, 1] != s[, 1]) > 0)
            {
              next;
            }
            gene_RPKM = cbind(gene_RPKM, s[, 2, drop = FALSE]);
            exp_names_gene = c(exp_names_gene, sample_TCGA_id);
          }
          
          time = proc.time() - time1;
          writeLines(paste("Downloaded - ", SpecificName, " - file ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        
        writeLines("Save data to local disk.");
        ID = str_locate_all(traverseResultFile, "_")[[1]];
        ID = ID[dim(ID)[1], 2];
        filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__gene.quantification__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
        write.table(rbind(c("Hybridization REF", exp_names_gene), c("Composite Element REF", column_gene), cbind(gene_left_column, gene_RPKM)), 
                    file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
        
        # For returning downloaded data
        dataIndex = dataIndex + 1;          
        downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names_gene), c("Composite Element REF", column_gene), cbind(gene_left_column, gene_RPKM));
        names(downloadedData)[dataIndex] = paste(SpecificName, "gene.quantification", sep = "__");
        
      }
    }
  }
  
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
  
  # Return downloaded data
  if (length(downloadedData) > 0)
  {
    for (i in 1:length(downloadedData))
    {
      rownames(downloadedData[[i]]) = NULL;
      colnames(downloadedData[[i]]) = NULL;
    }
  }
  downloadedData;
}




DownloadmiRNASeqData <- function(traverseResultFile, saveFolderName, cancerType, assayPlatform = "miRNASeq", tissueType = NULL, inputPatientIDs = NULL, outputFileName = "")
{
  options(warn=-1);
  
  writeLines("**********************************************************************************");
  writeLines("");
  writeLines(paste("Download miRNA-seq data of ", cancerType, " patients.", sep = ""));
  
  # Check whether specific TCGA patient IDs are inputted. 
  if (!is.null(inputPatientIDs))
  {
    inputPatientIDs = toupper(gsub(pattern = "\\.", replacement = "-", x = inputPatientIDs));
  }  
  
  if ((outputFileName != "") & (!is.null(outputFileName)))
  {
    outputFileName = paste(outputFileName, "__", sep = "");
  }
  
  # For returning downloaded data
  downloadedData = vector("list", 0);     
  dataIndex = 0;     
  
  # load directory traverse result and create folder to store output files
  writeLines("Load information of TCGA data files.");
  file_url=c()
  upper_file_url=c()
  load(traverseResultFile);
  if (assayPlatform == "miRNASeq")
  {
    platform = c("illuminaga_mirnaseq", "illuminahiseq_mirnaseq");
  }
  dir.create(path = saveFolderName, recursive = TRUE);
  miRNALevel3ID = grep(pattern = toupper("miRNASeq\\.Level_3"), x = upper_file_url, ignore.case = FALSE);
  
  for (IDpl in 1:length(platform))
  {
    SpecificName = paste(cancerType, "__", "bcgsc.ca", "__", platform[IDpl], sep = "");
    SpecificID = grep(pattern = toupper(paste("/", cancerType, "/cgcc/bcgsc\\.ca/", platform[IDpl], sep = "")), x = upper_file_url, ignore.case = FALSE);
    ind = SpecificID[grepEnd(pattern = toupper("sdrf\\.txt"), x = upper_file_url[SpecificID], ignore.case = FALSE)];
    if (length(ind) == 0)
    {
      next;
    }
    if (length(ind) > 1)
    {
      URL = GetNewestURL(AllURL = file_url[ind]);
    }else{
      URL = file_url[ind];
    }  
    downloadResult = urlReadTable(url = URL);
    if (downloadResult$errorFlag != 0)
    {
      next;
    }
    sdrf = toupper(downloadResult$data);  
    
    # Process SDRF file, identify the columns of level 3 data file name and TCGA sample barcodes.
    level_3_filename_column = max(grep(pattern = "Derived Data File", x = sdrf[1, ], ignore.case = TRUE));
    DataLevelColID = max(grep(pattern = "Comment \\[TCGA Data Level]", x = sdrf[1, ], ignore.case = TRUE));
    ExtractNameColID = min(grep(pattern = "Comment \\[TCGA Barcode]", x = sdrf[1, ], ignore.case = TRUE));
    RefGenomeColID = grep(pattern = "Comment \\[Genome reference]", x = sdrf[1, ], ignore.case = TRUE);
    if (length(ExtractNameColID) == 0)
    {
      ExtractNameColID = min(grep(pattern = "Extract Name", x = sdrf[1, ], ignore.case = TRUE));
    }
    colnames(sdrf) = sdrf[1, ];
    sdrf = unique(sdrf[2:dim(sdrf)[1], , drop = FALSE]);
    sdrf = sdrf[!duplicated(sdrf[, level_3_filename_column]), , drop = FALSE];
    
    Level3_ID = sort(union(which(sdrf[, DataLevelColID] == "LEVEL_3"), which(sdrf[, DataLevelColID] == "LEVEL 3")), decreasing = FALSE);
    Level3_ID = sort(intersect(Level3_ID, grep(pattern = toupper("mirna\\.quantification"), 
                                               x = sdrf[, level_3_filename_column], ignore.case = FALSE)), decreasing = FALSE);
    if (length(Level3_ID) == 0)
    {
      next;
    }
    sdrf = sdrf[Level3_ID, c(ExtractNameColID, level_3_filename_column, RefGenomeColID), drop = FALSE];    
    
    # If specific patient TCGA barcodes are inputted, only download the specified samples.
    if (!is.null(inputPatientIDs))
    {
      indInputPatientID = c();
      for (i in 1:length(inputPatientIDs))
      {
        indInputPatientID = c(indInputPatientID, grepBeginning(pattern = toupper(inputPatientIDs[i]), x = sdrf[, 1], ignore.case = FALSE));
      }
      if (length(indInputPatientID) == 0)
      {
        next;      
      }else{
        sdrf = sdrf[indInputPatientID, , drop = FALSE];
      }
    }
    
    # Download data of specified tissue
    if (!is.null(tissueType))
    {
      SampleType = cbind(Options = c("TP", "TR", "TB", "TRBM", "TAP", "TM", "TAM", "THOC", "TBM", "NB", "NT", "NBC", "NEBV", "NBM"),
                         Code = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14"));  
      sdrf = sdrf[substr(sdrf[, 1], 14, 15) %in% SampleType[SampleType[, "Options"] %in% tissueType, "Code"], , drop = FALSE];
    }
    
    if (dim(sdrf)[1] == 0)
    {
      next;
    }
    
    sdrfO = sdrf;
    for (RefGId in 1:2)
    {
      if (RefGId == 1)
      {
        sdrf = sdrfO[grep(pattern = toupper("NCBI36"), x = sdrfO[, 3], ignore.case = FALSE), , drop = FALSE];  
      }
      if (RefGId == 2)
      {
        sdrf = sdrfO[grep(pattern = toupper("GRCh37"), x = sdrfO[, 3], ignore.case = FALSE), , drop = FALSE];  
      }
      if (dim(sdrf)[1] == 0)
      {
        next;
      }
      
      exp_names = NULL;
      gene_left_column = NULL;
      gene_RPM = NULL;
      column_gene = NULL;
      for (i in 1:dim(sdrf)[1])
      {
        time1 = proc.time();
        
        sample_TCGA_id = sdrf[i, 1];
        ind = intersect(miRNALevel3ID, SpecificID[grepEnd(pattern = sdrf[i, 2], x = upper_file_url[SpecificID], ignore.case = FALSE)]);
        if (length(ind) == 0)
        {
          next;
        }
        if (length(ind) > 1)
        {
          URL = GetNewestURL(AllURL = file_url[ind]);
        }else{
          URL = file_url[ind];
        }  
        downloadResult = urlReadTable(url = URL);
        if (downloadResult$errorFlag != 0)
        {
          next;
        }
        s = downloadResult$data[2:dim(downloadResult$data)[1], , drop = FALSE];
        I_order_probes = order(s[, 1], decreasing = FALSE);
        s = s[I_order_probes, , drop = FALSE];
        if (is.null(gene_left_column))
        {
          gene_left_column = s[, 1, drop = FALSE];
          exp_names = c(sample_TCGA_id, sample_TCGA_id);
          gene_RPM = s[, 2:3, drop = FALSE];
          column_gene = c("read_count", "reads_per_million_miRNA_mapped");
        }else{
          if (sum(gene_left_column[, 1] != s[ ,1]) > 0)
          {
            next;
          }
          exp_names = c(exp_names, sample_TCGA_id, sample_TCGA_id);
          gene_RPM = cbind(gene_RPM, s[, 2:3, drop = FALSE]);
          column_gene = c(column_gene, "read_count", "reads_per_million_miRNA_mapped");
        }
        
        time = proc.time() - time1;
        if (RefGId == 1)
        {
          writeLines(paste("Downloaded - ", SpecificName, " - NCBI36 - sample ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        if (RefGId == 2)
        {
          writeLines(paste("Downloaded - ", SpecificName, " - GRCh37 - sample ", i, " out of ", dim(sdrf)[1], ". ", round(time[3], digits = 1), " seconds elapsed.", sep = ""));
        }
        
      }
      
      if (!is.null(gene_RPM))
      {
        writeLines("Save data to local disk.");
        ID = str_locate_all(traverseResultFile, "_")[[1]];
        ID = ID[dim(ID)[1], 2];
        if (RefGId == 1)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__NCBI36__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names), c("miRNA_ID", column_gene), cbind(gene_left_column, gene_RPM)), file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex+1;     
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names), c("miRNA_ID", column_gene), cbind(gene_left_column, gene_RPM));
          rownames(downloadedData[[dataIndex]]) = NULL;  
          colnames(downloadedData[[dataIndex]]) = NULL;          
          names(downloadedData)[dataIndex] = paste(SpecificName, "__NCBI36", sep = "");
        }
        if (RefGId == 2)
        {
          filename = paste(saveFolderName, "/", outputFileName, SpecificName, "__GRCh37__", substr(traverseResultFile, ID+1, nchar(traverseResultFile)-4), ".txt", sep = "");
          write.table(rbind(c("Hybridization REF", exp_names), c("miRNA_ID", column_gene), cbind(gene_left_column, gene_RPM)), file = filename, quote = FALSE, sep = "\t", na = "", col.names = FALSE, row.names = FALSE);  
          
          # For returning downloaded data
          dataIndex = dataIndex+1;     
          downloadedData[[dataIndex]] = rbind(c("Hybridization REF", exp_names), c("miRNA_ID", column_gene), cbind(gene_left_column, gene_RPM));
          rownames(downloadedData[[dataIndex]]) = NULL;  
          colnames(downloadedData[[dataIndex]]) = NULL;
          names(downloadedData)[dataIndex] = paste(SpecificName, "__GRCh37", sep = "");
        }
      }
      
    }
  }
  
  writeLines("");
  writeLines("**********************************************************************************");
  writeLines("\n");
  
  options(warn=0);
  
  downloadedData;
}



######################### Auxiliary Functions of Module A #####################################################

# urlReadTable is the function to read a data table from a website and pass it to a variable.

# Input arguments:
# url: URl of the website from which data table will be obtained.

# Output argument:
# data: a character matrix holding the data table obtained from website. 

urlReadTable <- function(url)
{
  
  data = try(content(GET(url), as = "text"), silent = TRUE);
  if (class(data) == "try-error")
  {
    return(list(data = data, errorFlag = 1));    
  }
  if (length(grep(pattern = "HTTP 404: Page Not Found\n\nThe page you requested was not found.", x = data, ignore.case = TRUE)) > 0)
  {
    return(list(data = data, errorFlag = 2)); 
  }
  data = gsub(pattern = "\r", replacement = "", x = data);
  data = as.matrix(read.table(text = data, sep="\t", fill = TRUE, quote = NULL, check.names = FALSE));
    
  return(list(data = data, errorFlag = 0));
}



# downloadFile is a function to download content from a website and save it as a local file.

# Input arguments:
# url: URl of the website whose content to be obtained.
# saveFileName: path and name of the file to store the web content.
##########delete downloadFile as we may not need it.



# grepEnd is a function similar to grep but identifies the strings with pattern at the end of the strings.
grepEnd <- function(pattern, x, ignore.case = FALSE)
{
  ind = grep(pattern = pattern, x = x, ignore.case = ignore.case);
  if (ignore.case)
  {
    ind = ind[nchar(x[ind]) == sapply(str_locate_all(toupper(x[ind]), pattern = toupper(pattern)), function(y){y[dim(y)[1], 2]})];    
  }else{
    ind = ind[nchar(x[ind]) == sapply(str_locate_all(x[ind], pattern = pattern), function(y){y[dim(y)[1], 2]})];    
  }
}



# grepBeginning is a function similar to grep but identifies the strings with pattern at the Beginning of the strings.
grepBeginning <- function(pattern, x, ignore.case = FALSE)
{
  ind = grep(pattern = pattern, x = x, ignore.case = ignore.case);
  if (ignore.case)
  {
    ind = ind[rep(1, length(ind)) == sapply(str_locate_all(toupper(x[ind]), pattern = toupper(pattern)), function(y)y[1, 1])];    
  }else{
    ind = ind[rep(1, length(ind)) == sapply(str_locate_all(x[ind], pattern = pattern), function(y)y[1, 1])];    
  }
}



# This function analyzes a TCGA webpage and identify all files and directories on it. 
# The URLs of files and directories are returned using two character vectors.

# Input arguments:
# TCGA_link: a string of URL for the TCGA webpage to be analyzed.

# Output arguments:
# file_url: a character vector, including the URLs of all files on the webpage.
# dir_url: a character vector, including the URLs of all directories on the webpage.


############Delete getTCGA_URL as it is only needed for TraverseAllDirectories to get the general directories. We are using the downloaded file instead.


# This function selects the newest file URL from all the input URLs by considering the
# last folder name tail numbers. The number format is *.*.*

GetNewestURL <- function(AllURL)
{
  SeriesNum = matrix(rep("", length(AllURL)*3), length(AllURL), 3);
  NumLength = rep(0, 3);
  SN = rep("", length(AllURL));
  for (i in 1:length(AllURL))
  {
    Str = AllURL[i];
    SepID = str_locate_all(string = Str, pattern = "/")[[1]];
    SepID = SepID[(dim(SepID)[1]-1):dim(SepID)[1], 1];
    Str = substr(Str, SepID[1]+1, SepID[2]-1);
    Str = strsplit(Str, split = "\\.")[[1]];
    SeriesNum[i, ] = Str[(length(Str)-2) : length(Str)];
    NumLength[1] = max(NumLength[1], nchar(SeriesNum[i, 1]));
    NumLength[2] = max(NumLength[2], nchar(SeriesNum[i, 2]));
    NumLength[3] = max(NumLength[3], nchar(SeriesNum[i, 3]));    
  }
  
  for (i in 1:length(AllURL))
  {
    for (j in 1:3)
    {
      if (nchar(SeriesNum[i, j]) < NumLength[j])
      {
        SeriesNum[i, j] = paste(paste(rep("0", NumLength[j] - nchar(SeriesNum[i, j])), collapse = ""), SeriesNum[i, j], sep = "");
      }
    }
    SN[i] = paste(SeriesNum[i, ], collapse = "");
  }
  AllURL[which.max(as.numeric(SN))];
}




#######################################################end TCGAAssembler - ModuleA.R######################################
#######################################################TCGAAssembler - MOduleB.R########################################
######################### Main Functions of Module B #############################################################

CombineMultiPlatformData<-function(inputDataList, combineStyle = "Intersect")
{
  options(warn=-1);
  
  for (i in 1:length(inputDataList))
  {
    
    # Keep only one sample for a tissue type of a patient. Usually, there is only one sample of a tissue type of a patient existing in data.
    inputDataList[[i]]$Data = ToPatientData(inputDataList[[i]]$Data);
    if (i == 1)
    {
      sample = colnames(inputDataList[[i]]$Data);
    }
    # For each genomic feature, keep only one row of data.
    Result = CombineRedundantFeature(Data = inputDataList[[i]]$Data, Des = inputDataList[[i]]$Des);
    
    inputDataList[[i]]$Des = Result$Des;
    inputDataList[[i]]$Data = Result$Data;
    if (combineStyle == "Intersect")
    {
      sample = sort(intersect(sample, colnames(inputDataList[[i]]$Data)));
    }
    if (combineStyle == "Union")
    {
      sample = sort(union(sample, colnames(inputDataList[[i]]$Data)));
    }    
    if (dim(inputDataList[[i]]$Des)[2] == 1)
    {
      inputDataList[[i]]$Des = cbind(inputDataList[[i]]$Des, Description = cbind(rep("", dim(inputDataList[[i]]$Des)[1])));
    }
    else
    {
      if ((dim(inputDataList[[i]]$Des)[2] == 3) && (inputDataList[[i]]$dataType == "CNA"))
      {
        inputDataList[[i]]$Des = cbind(inputDataList[[i]]$Des[, 1, drop = FALSE], 
                                       Description = paste(inputDataList[[i]]$Des[, 2], inputDataList[[i]]$Des[, 3], sep = ""));
      }
    }
    
    if (inputDataList[[i]]$dataType == "miRNAExp")
    {
      for (kk in 1:dim(inputDataList[[i]]$Des)[1])
      {
        inputDataList[[i]]$Des[kk, 1] = paste(substr(inputDataList[[i]]$Des[kk, 1], 5, 7), substr(inputDataList[[i]]$Des[kk, 1], 9, 100), sep = "");
      }
      inputDataList[[i]]$Des[, 1] = toupper(inputDataList[[i]]$Des[, 1]);
    }
    
    # for combining methylation data at CpG site level.
    if (inputDataList[[i]]$dataType == "Methylation")
    {    
      if (sum(toupper(c("REF", "GeneSymbol", "ChromosomeID", "CoordinateID")) %in% toupper(colnames(inputDataList[[i]]$Des))) == 4)
      {
        inputDataList[[i]]$Des = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
                                       Description = paste(inputDataList[[i]]$Des[, "REF"], inputDataList[[i]]$Des[, "ChromosomeID"], 
                                                           inputDataList[[i]]$Des[, "CoordinateID"], sep = "|"));
        inputDataList[[i]]$Des[, "Description"] = gsub(" ", "", inputDataList[[i]]$Des[, "Description"]);
      }
    }
    
    inputDataList[[i]]$Des = switch(inputDataList[[i]]$dataType,
                                    GeneExp = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
                                                    Platform = rep("GE", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]),
                                    ProteinExp = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
                                                       Platform = rep("PE", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]),                                    
                                    Methylation = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
                                                        Platform = rep("ME", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]),                                    
                                    CNA = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE], 
                                                Platform = rep("CN", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]), 
                                    miRNAExp = cbind(inputDataList[[i]]$Des[, "GeneSymbol", drop = FALSE],
                                                     Platform = rep("miRExp", dim(inputDataList[[i]]$Des)[1]), Description = inputDataList[[i]]$Des[, 2]));
    
    # check NA gene. which(inputDataList[[i]]$Des[, "GeneSymbol"] != "NA")
    ID = intersect(which(!is.na(inputDataList[[i]]$Des[, "GeneSymbol"])), 
                   which(inputDataList[[i]]$Des[, "GeneSymbol"] != "?"));
    inputDataList[[i]]$Des = inputDataList[[i]]$Des[ID, , drop = FALSE];
    inputDataList[[i]]$Data = inputDataList[[i]]$Data[ID, , drop = FALSE];    
    
  }
  
  for (i in 1:length(inputDataList))
  {
    # Get the data of samples that should be kept
    if (combineStyle == "Intersect")
    {
      inputDataList[[i]]$Data = inputDataList[[i]]$Data[, sample, drop = FALSE];  
    }
    if (combineStyle == "Union")
    {
      tempData = matrix(NA, dim(inputDataList[[i]]$Data)[1], length(sample));
      colnames(tempData) = sample;
      tempData[, colnames(inputDataList[[i]]$Data)] = inputDataList[[i]]$Data;
      inputDataList[[i]]$Data = tempData;  
    }
    rownames(inputDataList[[i]]$Data) = NULL;
    rownames(inputDataList[[i]]$Des) = NULL;  
  }
  
  # Combine the datasets into matrix format
  Data = inputDataList[[1]]$Data;
  Des = inputDataList[[1]]$Des;  
  for (i in 2:length(inputDataList))
  {
    Data = rbind(Data, inputDataList[[i]]$Data);
    Des = rbind(Des, inputDataList[[i]]$Des);    
  } 
  
  OrderID = order(as.character(Des[, "GeneSymbol"]), as.character(Des[, "Platform"]), 
                  as.character(Des[, "Description"]), na.last = TRUE, decreasing = FALSE);
  Data = Data[OrderID, , drop = FALSE];
  Des = Des[OrderID, , drop = FALSE];
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  Result = list(Des = Des, Data = Data);
  
  options(warn=0);
  
  Result;
}

##########delete ExtractTissueSpecificSamples


ProcessRNASeqData<-function(inputFilePath, outputFileName, outputFileFolder, dataType, verType)
{
  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  # Read in data.
  InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.strings = "", stringsAsFactors = FALSE, quote = "", check.names = FALSE);
  InData = InData[2:dim(InData)[1], , drop = FALSE];  
  if ((dataType == "GeneExp") & (verType == "RNASeqV1"))
  {
    REF = InData[, 1];
    RPKM = as.matrix(InData[, seq(4, dim(InData)[2], 3), drop = FALSE]);
    mode(RPKM) = "numeric";
    GeneSymbol = sapply(strsplit(REF, split = "\\|"), function(x)x[1]);
    EntrezID = sapply(strsplit(REF, split = "\\|"), function(x)x[2]);
    Des = cbind(GeneSymbol = GeneSymbol, EntrezID = EntrezID);
    Data = RPKM;
    #Check and Correct gene symbol
    Des = CheckGeneSymbol(Des);    
  }
  if (dataType == "ExonExp")
  {
    REF = InData[, 1];
    RPKM = as.matrix(InData[, seq(4, dim(InData)[2], 3), drop = FALSE]);
    mode(RPKM) = "numeric";
    Des = cbind(ExonID = REF);
    Data = RPKM;    
  }
  if ((dataType == "GeneExp") & (verType == "RNASeqV2"))
  {
    REF = InData[, 1];
    NormalizedCount = as.matrix(InData[, 2:dim(InData)[2], drop = FALSE]);
    mode(NormalizedCount) = "numeric";
    GeneSymbol = sapply(strsplit(REF, split = "\\|"), function(x)x[1]);
    EntrezID = sapply(strsplit(REF, split = "\\|"), function(x)x[2]);
    Des = cbind(GeneSymbol = GeneSymbol, EntrezID = EntrezID);
    Data = NormalizedCount;    
    #Check and Correct gene symbol
    Des = CheckGeneSymbol(Des);     
  }
  if (verType == "Microarray")
  {
    dataType = "GeneExp";
    REF = InData[, 1];
    Data = as.matrix(InData[, 2:dim(InData)[2], drop = FALSE]);
    mode(Data) = "numeric";
    Des = cbind(GeneSymbol = REF, EntrezID = rep("", length(REF)));
    #Check and Correct gene symbol
    Des = CheckGeneSymbol(Des);     
  }
  
  if (dataType == "GeneExp")
  {
    # Draw and save a box plot
    png(filename = paste(outputFileFolder, "/", outputFileName, "__boxplot.png", sep = ""), width = 30*dim(Data)[2]+300, height = 1500, units = "px");
    par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
    nonNaID = which(!is.na(Data));
    if ((sum(Data[nonNaID]<0)==0) && (max(Data[nonNaID])>50))
    {
      boxplot(log2(Data), main = "Boxplot drawn based on log2 tranformed data.");
    } else {
      boxplot(Data);
    }
    dev.off(); 
  }
  
  # save data files
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, ".rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, ".txt", sep = ""), quote = FALSE, 
              sep = "\t", na = "", col.names = TRUE, row.names = FALSE); 
  
  options(warn=0);
  
  return(list(Data = Data, Des = Des)); 
}



ProcessmiRNASeqData<-function(inputFilePath, outputFileName, outputFileFolder, fileSource = "TCGA-Assembler")
{
  options(warn=-1);
  
  dir.create(path = outputFileFolder, recursive = TRUE);
  
  # Read in data.
  InData = read.table(file = inputFilePath, header = TRUE, sep = "\t", na.strings = "", quote = "", stringsAsFactors = FALSE, check.names = FALSE);
  InData = InData[2:dim(InData)[1], , drop = FALSE];
  
  # divide read cout data and RPM data
  REF = InData[, 1];  
  if (fileSource == "Firehose")
  {
    Count = as.matrix(InData[, seq(2, dim(InData)[2], 3), drop = FALSE]);
    RPM = as.matrix(InData[, seq(3, dim(InData)[2], 3), drop = FALSE]);
  }
  if (fileSource == "TCGA-Assembler")
  {
    Count = as.matrix(InData[, seq(2, dim(InData)[2], 2), drop = FALSE]);
    RPM = as.matrix(InData[, seq(3, dim(InData)[2], 2), drop = FALSE]);    
  }
  #  rownames(Count) = REF;  
  mode(Count) = "numeric";
  #  rownames(RPM) = REF;
  mode(RPM) = "numeric";
  
  #   #   Draw and save box plot
  #   png(filename = paste(outputFileFolder, "/", outputFileName, "__ReadCount.boxplot.png", sep = ""), width = 30*dim(Count)[2]+300, height = 1500, units = "px");
  #   par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  #   boxplot(log2(Count));
  #   dev.off();
  #   png(filename = paste(outputFileFolder, "/", outputFileName, "__RPM.boxplot.png", sep = ""), width = 30*dim(RPM)[2]+300, height = 1500, units = "px");
  #   par(las = 2, mai = c(4.5, 1, 0.5, 1), cex = 1.5);
  #   boxplot(log2(RPM));
  #   dev.off();  
  
  #save data files
  Des = cbind(GeneSymbol = REF);
  Data = Count;
  rownames(Data) = NULL;
  rownames(Des) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, "__ReadCount.rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, "__ReadCount.txt", sep = ""), quote = FALSE, 
              sep = "\t", na = "", col.names = TRUE, row.names = FALSE); 
  Data = RPM; 
  rownames(Data) = NULL;
  save(Des, Data, file = paste(outputFileFolder, "/", outputFileName, "__RPM.rda", sep = ""));
  write.table(cbind(Des, Data), file = paste(outputFileFolder, "/", outputFileName, "__RPM.txt", sep = ""), quote = FALSE, 
              sep = "\t", na = "", col.names = TRUE, row.names = FALSE); 
  
  options(warn=0);
  
  return(list(Data = Data, Des = Des)); 
}



CheckGeneSymbol<-function(Des)
{
  data("hgnc.table", package="HGNChelper", envir=environment());
  hgnc.table = rbind(hgnc.table, c("13-SEP", "SEPT7P2"));
  rID = intersect(which(toupper(hgnc.table[, "Symbol"]) == toupper("NCRNA00185")), 
                  which(toupper(hgnc.table[, "Approved.Symbol"]) == toupper("TTTY14")));
  hgnc.table = hgnc.table[sort(setdiff(1:dim(hgnc.table)[1], rID)), , drop = FALSE];
  hgnc.table = hgnc.table[which(!is.na(hgnc.table[, "Approved.Symbol"])), , drop = FALSE];
  regex = "[0-9]\\-(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)|[0-9]\\.[0-9][0-9]E\\+[[0-9][0-9]";
  MonthID = grep(pattern = regex, hgnc.table[, 1], ignore.case = TRUE);
  MonthMappingTable = hgnc.table[MonthID, , drop = FALSE];
  Des[, "GeneSymbol"] = sub("(.*C[0-9XY]+)ORF(.+)", "\\1orf\\2", Des[, "GeneSymbol"]);
  ID = intersect(which(!Des[, "GeneSymbol"] %in% hgnc.table[, "Approved.Symbol"]), which(Des[, "GeneSymbol"] %in% hgnc.table[, "Symbol"]));
  DesOrg = Des;
  if (length(ID) > 0)
  {
    Des[ID, "GeneSymbol"] = sapply(Des[ID, "GeneSymbol"], function(x)paste(hgnc.table[hgnc.table[, "Symbol"] == x, "Approved.Symbol"], collapse = "___"));
  }
  #  writeLines("Changed genes are");
  #  print(cbind(DesOrg[ID, "GeneSymbol"], Des[ID, "GeneSymbol"]));  
  ID = intersect(which(!Des[, "GeneSymbol"] %in% hgnc.table[, "Approved.Symbol"]), which(toupper(Des[, "GeneSymbol"]) %in% toupper(MonthMappingTable[, "Symbol"])));
  if (length(ID) > 0)
  {
    Des[ID, "GeneSymbol"] = sapply(Des[ID, "GeneSymbol"], function(x)paste(MonthMappingTable[toupper(MonthMappingTable[, "Symbol"]) == toupper(x), "Approved.Symbol"], collapse = "___"));
  }
  #  writeLines("Changed genes are");
  #  print(cbind(DesOrg[ID, "GeneSymbol"], Des[ID, "GeneSymbol"]));  
  Des;
}

#delete ProcessRPPADataWithGeneAnnotation

######################### Auxiliary Functions of Module B #############################################################

ToPatientData<-function(TCGAData)
{
  TCGABarcode = colnames(TCGAData);
  TCGABarcode = sapply(strsplit(TCGABarcode, split = "-"), function(x){
    Str = paste(x[1], x[2], x[3], substr(x[4], 1, 2), sep = "-");
    Str;
  });
  DuplicatedLabel = duplicated(TCGABarcode);
  ID = which(DuplicatedLabel == FALSE);
  TCGAData = TCGAData[, ID, drop = FALSE];
  TCGABarcode = TCGABarcode[ID];
  colnames(TCGAData) = TCGABarcode;
  TCGAData;
}



CombineRedundantFeature <- function(Data, Des)
{
  RowName = Des[, 1];
  if (dim(Des)[2] > 1)
  {
    for (i in 2:dim(Des)[2])
    {
      RowName = paste(RowName, Des[, i], sep = "||");
    }
  }
  
  UniqueRowName = unique(RowName);
  if (length(RowName) == length(UniqueRowName))
  {
    Result = list(Data = Data, Des = Des);
  }
  else
  {
    IDKeep = c();
    for (i in 1:length(UniqueRowName))
    {
      IDi = which(RowName == UniqueRowName[i]);
      if (length(IDi) > 1)
      {
        Data[IDi[1], ] = apply(Data[IDi, , drop = FALSE], 2, mean, na.rm = TRUE);
      }
      IDKeep = c(IDKeep, IDi[1]);
    }
    IDKeep = sort(IDKeep);
    Result = list(Data = Data[IDKeep, , drop = FALSE], Des = Des[IDKeep, , drop = FALSE]);
  }
  Result;
}


###delete normalize.quantiles


#####################
############Delete VCwebContent


#######################################################end TCGAAssembler ModulB.R###########################




#' Get data from TCGA
#' 
#' Download matched miRNA and mRNA expression profiles from TCGA using TCGA Assembler. In order to this function to work, users need to download the
#' the file "DirectoryTraverseResult_Jul-08-2014.rda" from nugget.unisa.edu.au/Thuc/miRLAB/ and place it in the working directory 

#' @param CancerName The name of cancer using abbreviation of cancer name as one of the followings:
#' \itemize{
#' \item Adrenocortical carcinoma: ACC
#' \item Bladder urothelial carcinoma: BLCA
#' \item Breast invasive carcinoma: BRCA
#' \item Cervical squamous cell carcinoma and endocervical adenocarcinoma: CESC
#' \item Colon adenocarcinoma: COAD
#' \item Lymphoid neoplasm diffuse large B-cell lymphoma: DLBC
#' \item Esophageal carcinoma: ESCA
#' \item Glioblastoma multiforme: GBM
#' \item Head and neck squamous cell carcinoma: HNSC
#' \item Kidney chromophobe: KICH
#' \item Kidney renal clear cell carcinoma: KIRC
#' \item Kidney renal papillary cell carcinoma: KIRP
#' \item Acute myeloid leukemia: LAML
#' \item Brain lower grade glioma: LGG
#' \item Liver hepatocellular carcinoma: LIHC
#' \item Lung adenocarcinoma: LUAD
#' \item Lung squamous cell carcinoma: LUSC
#' \item Ovarian serous cystadenocarcinoma: OV
#' \item Pancreatic adenocarcinoma: PAAD
#' \item Prostate adenocarcinoma: PRAD
#' \item Rectum adenocarcinoma: READ
#' \item Sarcoma: SARC
#' \item Skin cutaneous melanoma: SKCM
#' \item Stomach adenocarcinoma: STAD
#' \item Thyroid carcinoma: THCA
#' \item Uterine corpus endometrial carcinoma: UCEC
#' \item Uterine carcinosarcoma: UCS
#' }
#' @return Dataset of matched miRNA and mRNA expression profiles 
#' @examples 
#' \dontrun{
#' BreastCancer=getData("BRCA") 
#' }
#' @references 
#' Zhu, Y., Qiu, P., and Ji, Y. (2014). TCGA-Assembler: open-source software for retrieving and processing TCGA data. Nat. Methods, 11, 599-600.
#' @export
getData<-function(CancerName){
#rm(list = ls()); # Clear workspace
#source("Module_A.r"); # Load Module A functions
#source("Module_B.r"); # Load Module B functions

# Get the URLs of public data files of all cancer types and all assay platforms.
# TraverseAllDirectories(entryPoint = "https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/", fileLabel = "DirectoryTraverseResult");

# Download and process miRNA-seq data of all CancerName patient samples.
miRNASeqRawData = DownloadmiRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./UserData/RawData.TCGA-Assembler", cancerType = CancerName, assayPlatform = "miRNASeq");
miRNADataName=substr(colnames(data.frame(miRNASeqRawData[1]))[1],1,nchar(colnames(data.frame(miRNASeqRawData[1]))[1])-2);
InputmiRNAfilename = paste("./UserData/RawData.TCGA-Assembler/",miRNADataName,"__Jul-08-2014",".txt",sep="");
OutputmiRNAfilename = paste(CancerName,"_illuminahiseq_mirnaseq",sep="");
miRNASeqData = ProcessmiRNASeqData(inputFilePath = InputmiRNAfilename, outputFileName = OutputmiRNAfilename, outputFileFolder = "./UserData/ProcessedData.TCGA-Assembler");

# Download and process normalized gene expression data of all CancerName patient samples
RNASeqRawData = DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./UserData/RawData.TCGA-Assembler", cancerType = CancerName, assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results");
RNADataName=substr(colnames(data.frame(RNASeqRawData[1]))[1],1,nchar(colnames(data.frame(RNASeqRawData[1]))[1])-2);
InputRNAfilename = paste("./UserData/RawData.TCGA-Assembler/",RNADataName,"__Jul-08-2014",".txt",sep="");
OutputRNAfilename = paste(CancerName,"_illuminahiseq_rnaseqv2__GeneExp",sep="");
GeneExpData = ProcessRNASeqData(inputFilePath = InputRNAfilename, outputFileName = OutputRNAfilename, outputFileFolder = "./UserData/ProcessedData.TCGA-Assembler", dataType = "GeneExp", verType = "RNASeqV2");

# Put multi-modal data in a vector of list objects to be inputted into CombineMultiPlatformData function.
inputDataList = vector("list", 2);
inputDataList[[1]] = list(Des = miRNASeqData$Des, Data = miRNASeqData$Data, dataType = "miRNAExp");
inputDataList[[2]] = list(Des = GeneExpData$Des, Data = GeneExpData$Data, dataType = "GeneExp");

# Merge multi-platform data using Intersect approach.
MergedData = CombineMultiPlatformData(inputDataList = inputDataList);
# Merge multi-platform data using Union approach.
# MergedData = CombineMultiPlatformData(inputDataList = inputDataList, combineStyle = "Union");
   
return(MergedData)

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
#' #' @param r The rate threshold to filter the records (genes). Genes with more than r\% missing data will be removed.
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
#' @param topkmiR the number of miRNAs that we would like to extract, e.g. top 50 miRNAs.
#' @param topkmR the number of mRNAs that we would like to extract, e.g. top 2000 mRNAs.
#' @return the dataset that includes differentially expressed miRNAs and mRNAs. columns are miRNAs and mRNAs  and rows are samples
#' @references
#' Smyth, G.K. (2005). Limma: linear models for microarray data. In Bioinformatics and computational biology solutions using R and Bioconductor (pp. 397-420). Springer New York.
#' @export
DiffExpAnalysis=function(miR1, miR2, mR1, mR2,topkmiR, topkmR){
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
	miRresults=topTable(miRfit2, number= topkmiR, sort.by="p", adjust.method="BH")
	write.csv(miRresults, file="DiffExpmiR.csv")
	########## mR ############
	mRfit=lmFit(mR, design)
	contrast.matrix=makeContrasts(NormalvCancer=Normal - Cancer, levels=design)
	mRfit2=contrasts.fit(mRfit, contrast.matrix)
	mRfit2=eBayes(mRfit2)
	mRresults=topTable(mRfit2, number= topkmR, sort.by="p", adjust.method="BH")
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

				