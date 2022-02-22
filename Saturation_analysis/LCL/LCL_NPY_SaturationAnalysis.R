#This is a script to determine whether I have saturated the analysis. 
#Start with LCLs X-responsive genes
inputValue <- commandArgs(trailingOnly=TRUE)
# inputValue <- c(59,1)
print(as.character(inputValue[1]))
print(as.character(inputValue[2]))

set.seed(1)

myPath <- #PATH TO GITHUB FILES 

### 1. Bring in sample tables for aneuploidy samples ###
suppressPackageStartupMessages(library("limma"))

#read in the metadata table
metadata <-read.delim(paste0(myPath,"/Linear_regressions/metadata.txt"), stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample
#subtract 1 from X count to model number of additional X chromosomes
metadata$x_count <- metadata$x_count - 1
#restrict to male samples only for NPY analysis
metadata_lcl_male <- metadata[metadata$karyotype %in% c("XXXXY","XXXY","XXY","XXYY","XY","XYY","XYYYY") & metadata$cell_type == "LCL" & metadata$Technical_replicate == "N",]
rownames(metadata_lcl_male) <- metadata_lcl_male$sample
#subtract 1 from Y count to model number of additional Y chromosomes
metadata_lcl_male$y_count <- metadata_lcl_male$y_count - 1

#bring in normalized read counts table for expressed NPY genes
normCounts <- read.delim(file=paste0(myPath,"/Linear_regressions/LCLs/lcl_normCounts_expYgenes_males.txt"), check.names = FALSE)

#create results table
saturationResults <- NULL

#choose this many samples 100 times
i <- inputValue[1]

j <- 1
repeat {
  #randomly select samples
  mySamples <- sample(1:dim(metadata_lcl_male)[1], i, replace = FALSE)
  sample_table <- metadata_lcl_male[mySamples,]

  #Test whether matrix is full rank
  mod <- model.matrix(~batch_libprep + y_count + x_count, sample_table)
  if(is.fullrank(mod)){
    ##Do analysis of NPY genes##
    normCounts_expYgenes <- normCounts[, sample_table$sample]
    my_metadata <- data.frame(sample_table$x_count, sample_table$y_count, sample_table$batch_libprep)
    rownames(my_metadata) <- sample_table$sample
    
    myResults_npy <- data.frame(delEy = numeric(), ypval = numeric())
    for(ygene in rownames(normCounts_expYgenes)){
      
      #pull out data
      myGeneData <- t(normCounts_expYgenes[ygene,])
      all_Data <- data.frame(my_metadata, "ygene_expression" = myGeneData)
      colnames(all_Data) <- c("x_count", "y_count", "batch_libprep", "ygene_expression")
      
      #do a linear regression
      myFormula <- formula(ygene_expression ~ x_count + y_count + batch_libprep)
      myFormula_int <- formula(ygene_expression ~ x_count + y_count)
      
      mylm <- lm(myFormula, data=all_Data)
      mylm_int <- lm(myFormula_int, data=all_Data)
      
      #pull out y_count beta, p-value and add to table
      myInt <- summary(mylm_int)$coefficients["(Intercept)",1]
      YBeta <- summary(mylm)$coefficients["y_count",1]
      YPval <- summary(mylm)$coefficients["y_count",4]
      delEy <- YBeta/myInt
      myResults_npy[ygene,] <- c(delEy, YPval)
    }
    
    #Sometimes I get a NA results for everything in some rounds, want to redo this round if that's the case, so don't advance j.
    if(is.na(myResults_npy[1,"ypval"])){
      print("NA values, drawing samples again.")
    } else{
      #correct for multiple hypothesis testing and set signficance threshold to padj < 0.05 for the Y chromosome
      bh_adj_npy <- p.adjust(myResults_npy$ypval, method="BH")
      myResults_npy <- cbind(myResults_npy, yadj_pval = bh_adj_npy)
      sigResults_npy <- myResults_npy[myResults_npy$yadj_pval < 0.05,]
      numSig_npy <- dim(sigResults_npy)[1]
      sigResults_npy_pos <- sigResults_npy[sigResults_npy$delEy > 0,]
      numSig_npy_pos <- dim(sigResults_npy_pos)[1]
      sigResults_npy_neg <- sigResults_npy[sigResults_npy$delEy < 0,]
      numSig_npy_neg <- dim(sigResults_npy_neg)[1]
      
      #Record the number of significant Y genes in a table
      sigGenes <- c(numSig_npy, numSig_npy_pos, numSig_npy_neg) 
      saturationResults <- rbind(saturationResults, c(j,sigGenes , summary(sample_table$karyotype)))
      names(saturationResults) <- c("round","npy_sig", "npy_sig_up","npy_sig_down", names(summary(sample_table$karyotype)))
      
      
      #print the results out
      print(paste0("Finished round ",as.character(i),",",as.character(j),". # of Y responsive genes: ", as.character(numSig_npy),"."))
      
      #advance j
      j <- j + 1
    }
  } else{
    print("Model not full rank, drawing samples again.")
  }
  
  if(j == as.numeric(inputValue[2]) + 1){
    break
  }
  
}
##outside for loop
saturationResults <- as.data.frame(saturationResults)
colnames(saturationResults) <- c("round","npy_sig", "npy_sig_up","npy_sig_down",
                                 names(summary(sample_table$karyotype)))

write.table(saturationResults, file=paste0("saturationResults_100iterations_male_", as.character(i),".txt"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep="\t")

sessionInfo()

