#This is a script to determine whether I have saturated the analysis.
#Start with LCLs X-responsive genes
inputValue <- commandArgs(trailingOnly = TRUE)
# inputValue <- c(45,1)
print(as.character(inputValue[1])) #number of samples
print(as.character(inputValue[2])) #number of iterations

#setpath
myPath <- #PATH TO GITHUB FILES 

set.seed(1)

## Load required functions 
suppressPackageStartupMessages(library(limma))

#read in the metadata table and restrict to the appropriate samples
metadata <- read.delim(paste0(myPath,"/Linear_regressions/metadata.txt"), stringsAsFactors = FALSE)
rownames(metadata) <- metadata$sample
metadata_tri21 <- metadata[metadata$karyotype %in% c("XX","XY","XY_t21","XX_t21") & metadata$cell_type == "LCL" & metadata$Technical_replicate  == "N",c("Chr21_count","sex","batch_libprep")]
metadata_tri21$Chr21_count <- metadata_tri21$Chr21_count - 2

#bring in normalized counts
normCounts <- read.delim(file=paste0(myPath,"/Linear_regressions/Chr21/lcl_normCounts_chr21_genes.txt"), check.names = FALSE)

#choose this many samples
i <- inputValue[1]
#create results table
saturationResults <- NULL

j <- 1
repeat {
  #randomly select samples
  mySamples <- sample(1:dim(metadata_tri21)[1], i, replace = FALSE)
  sample_table_tri21 <- metadata_tri21[mySamples,]
  
  #set up model matrix and check if full rank:
  mod_chr21 <- model.matrix(~ batch_libprep + Chr21_count + sex, sample_table_tri21)
  if (is.fullrank(mod_chr21)) {
    ### Linear modeling of chr21 genes
    normCounts_exp21genes <- normCounts[,rownames(sample_table_tri21)]
    my_metadata_tri21 <- data.frame(sample_table_tri21$Chr21_count, sample_table_tri21$sex, sample_table_tri21$batch_libprep)
    rownames(my_metadata_tri21) <- sample_table_tri21$sample
    
    myResults_chr21 <- data.frame(delE21 = numeric(), chr21pval = numeric())
    for(gene in rownames(normCounts_exp21genes)){
       #pull out data
      myGeneData <- t(normCounts_exp21genes[gene,])
      all_Data <- data.frame(my_metadata_tri21, "gene_expression" = myGeneData)
      colnames(all_Data) <- c("Chr21_count", "sex", "batch_libprep", "gene_expression")

      #do a linear regression
      myFormula <- formula(gene_expression ~ Chr21_count + sex + batch_libprep)
      myFormula_int <- formula(gene_expression ~ Chr21_count + sex)
      mylm <- lm(myFormula, data=all_Data)
      mylm_int <- lm(myFormula_int, data=all_Data)
      
      #pull out beta, p-value and add to table
      intercept <- (summary(mylm_int)$coefficients["(Intercept)",1])/2 #divide by two to get the value for one chr21
      chr21beta <- summary(mylm)$coefficients["Chr21_count",1]
      chr21pval <- summary(mylm)$coefficients["Chr21_count",4]
      delE21 <- chr21beta / intercept
      myResults_chr21[gene,] <- c(delE21, chr21pval)
    }
    if (is.na(myResults_chr21[1, "chr21pval"])) {
      print("NA values, drawing samples again.")
    } else{
      #correct for multiple hypothesis testing and set signficance threshold to padj < 0.05 for the X chromosome
      bh_adj_chr21 <-
        p.adjust(myResults_chr21$chr21pval, method = "BH")
      myResults_chr21 <-
        cbind(myResults_chr21, adj_pval = bh_adj_chr21)
      sigResults_chr21 <-
        myResults_chr21[myResults_chr21$adj_pval < 0.05, ]
      numSig_chr21 <- dim(sigResults_chr21)[1]
    }
    
    #Record the number of significant X or total sex chrom responsive genes in a table
    sigGenes <-
      c(numSig_chr21)
    saturationResults <-
      rbind(saturationResults, c(j, sigGenes , summary(sample_table_tri21$karyotype)))
    names(saturationResults) <-
      c("round",
        "chr21_sig",
        names(summary(sample_table_tri21$karyotype)))

    #print the results out
    print(
      paste0(
        "Finished round ",
        as.character(i),
        ",",
        as.character(j),
        ". # of Chr21 responsive genes: ",
        as.character(numSig_chr21)
      )
    )
    
    #advance j
    j <- j + 1
  }
  else{
    print("Model not full rank, drawing samples again.")
  }
  
  if (j == as.numeric(inputValue[2]) + 1) {
    break
  }
}


##outside for loop
saturationResults <- as.data.frame(saturationResults)
colnames(saturationResults) <- c("round", "chr21",names(summary(sample_table_tri21$karyotype)))
write.table(
  saturationResults,
  file = paste0("saturationResults_100iterations_", as.character(i), ".txt"),
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

sessionInfo()
    


