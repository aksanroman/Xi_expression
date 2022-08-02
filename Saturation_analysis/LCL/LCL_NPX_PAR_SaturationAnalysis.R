#This is a script to determine whether I have saturated the analysis.
#Start with LCLs X-responsive genes
inputValue <- commandArgs(trailingOnly = TRUE)
# inputValue <- c(106,1)
print(as.character(inputValue[1]))
print(as.character(inputValue[2]))

set.seed(1)

### 1. Bring in sample tables for aneuploidy samples ###
suppressPackageStartupMessages(library("limma"))

myPath <- #PATH TO GITHUB FILES 
  
PAR1_genes <- c("AKAP17A","ASMT","ASMTL","CD99","CRLF2","CSF2RA","DHRSX","GTPBP6","IL3RA","P2RY8","PLCXD1","PPP2R3B","SHOX","SLC25A6","ZBED1")
PAR2_genes <- c("VAMP7","WASH6P","IL9R","SPRY3")
XY_pairs <- c("PRKX","NLGN4X","TBL1X","AMELX","TMSB4X","TXLNG","EIF1AX","ZFX","USP9X","DDX3X","KDM6A","TSPYL2","KDM5C","RPS4X","RBMX","SOX3","HSFX1")

#read in the metadata table
metadata <-read.delim(paste0(myPath,"/Linear_regressions/metadata.txt"), stringsAsFactors = FALSE) 
rownames(metadata) <- metadata$sample
#subtract 1 from X count to model number of additional X chromosomes
metadata$x_count <- metadata$x_count - 1
#restrict to male samples only for NPY analysis
metadata_lcl <- metadata[metadata$karyotype %in% c("X","XX","XXX","XXXX","XXXXY","XXXY","XXY","XXYY","XY","XYY","XYYYY") & metadata$cell_type == "LCL" & metadata$Technical_replicate == "N",]

#bring in normalized read counts per sample/gene
normCounts <- read.delim(file=paste0(myPath,"/Linear_regressions/LCLs/LCL_normCounts_expXgenes.txt"), check.names = FALSE) 

#choose this many samples 100 times
i <- inputValue[1]

#create results table
saturationResults <- NULL

j <- 1
repeat {
  #randomly select samples
  mySamples <-sample(1:dim(metadata_lcl)[1], i, replace = FALSE)
  sample_table <- metadata_lcl[mySamples,]

  #Test whether matrix is full rank
  mod_npx <- model.matrix( ~ batch_libprep + y_count + x_count, sample_table)
  if (is.fullrank(mod_npx)) {
      ##Do analysis of NPX genes##
      normCounts_expXgenes <- normCounts[,sample_table$sample]
      my_metadata <- data.frame(sample_table$x_count, sample_table$y_count, sample_table$batch_libprep)
      rownames(my_metadata) <- sample_table$sample
      
      myResults_npx_par <- data.frame(delEx = numeric(), xpval = numeric(), delEy = numeric(), ypval=numeric())
      for (xgene in rownames(normCounts_expXgenes)) {
        #pull out data
        myGeneData <- t(normCounts_expXgenes[xgene, ])
        all_Data <- data.frame(my_metadata, "xgene_expression" = myGeneData)
        colnames(all_Data) <-
          c("x_count",
            "y_count",
            "batch_libprep",
            "xgene_expression")
        
        #do a linear regression
        myFormula <-
          formula(xgene_expression ~ x_count + y_count + batch_libprep)
        myFormula_int <- formula(xgene_expression ~ x_count + y_count )
        mylm <- lm(myFormula, data = all_Data)
        mylm_int <- lm(myFormula_int, data=all_Data)
        
        #pull out x_count beta, p-value and add to table
        XBeta <- summary(mylm)$coefficients["x_count", 1]
        XPval <- summary(mylm)$coefficients["x_count", 4]
        myInt <- summary(mylm_int)$coefficients["(Intercept)", 1]
        if(xgene == "XIST"){
          myInt <- myInt + XBeta
        }
        delEx <- XBeta / myInt
        YBeta <- summary(mylm)$coefficients["y_count", 1]
        YPval <- summary(mylm)$coefficients["y_count", 4]
        delEy <- YBeta / myInt
        myResults_npx_par[xgene, ] <- c(delEx, XPval, delEy, YPval)
      }

      if (is.na(myResults_npx_par[1, "xpval"])) {
        print("NA values, drawing samples again.")
      } else{
        #correct for multiple hypothesis testing and set signficance threshold to padj < 0.05 for the X chromosome
        bh_adj_npx_x <- p.adjust(myResults_npx_par$xpval, method = "BH")
        bh_adj_npx_y <- p.adjust(myResults_npx_par$ypval, method = "BH")
        myResults_npx_par <-
          cbind(myResults_npx_par, xadj_pval = bh_adj_npx_x, yadj_pval = bh_adj_npx_y)
        myResults_par <- myResults_npx_par[rownames(myResults_npx_par) %in% PAR1_genes,]
        myResults_npx <- myResults_npx_par[! rownames(myResults_npx_par) %in% c(PAR1_genes,PAR2_genes),]
        sigResults_npx <-
          myResults_npx[myResults_npx$xadj_pval < 0.05, ]
        numSig_npx <- dim(sigResults_npx)[1]
        sigResults_npx_pos <-
          sigResults_npx[sigResults_npx$delEx > 0, ]
        numSig_npx_pos <- dim(sigResults_npx_pos)[1]
        sigResults_npx_neg <-
          sigResults_npx[sigResults_npx$delEx < 0, ]
        numSig_npx_neg <- dim(sigResults_npx_neg)[1]
        sigResults_npx_pos_ypair <- sigResults_npx_pos[rownames(sigResults_npx_pos) %in% XY_pairs,]
        numSig_npx_pos_ypair <- dim(sigResults_npx_pos_ypair)[1]
        sigResults_npx_pos_noY <- sigResults_npx_pos[! rownames(sigResults_npx_pos) %in% XY_pairs,]
        numSig_npx_pos_noY <- dim(sigResults_npx_pos_noY)[1]
        
        sigResults_par <-
          myResults_par[myResults_par$xadj_pval < 0.05 & myResults_par$yadj_pval < 0.05, ]
        numSig_par <- dim(sigResults_par)[1]
        sigResults_par_pos <-
          sigResults_par[sigResults_par$delEx > 0, ]
        numSig_par_pos <- dim(sigResults_par_pos)[1]
        sigResults_par_neg <-
          sigResults_par[sigResults_par$delEx < 0, ]
        numSig_par_neg <- dim(sigResults_par_neg)[1]
        
        #Record the number of significant X or total sex chrom responsive genes in a table
        sigGenes <-
          c(
            numSig_npx,
            numSig_npx_pos,
            numSig_npx_pos_ypair,
            numSig_npx_pos_noY,
            numSig_npx_neg,
            numSig_par,
            numSig_par_pos,
            numSig_par_neg
          )
        saturationResults <-
          rbind(saturationResults, c(j, sigGenes , summary(sample_table$karyotype)))
        names(saturationResults) <-
          c(
            "round",
            "npx_sig",
            "npx_sig_up",
            "npx_sig_up_ypair",
            "npx_sig_up_noY",
            "npx_sig_down",
            "par_sig",
            "par_sig_up",
            "par_sig_down",
            names(summary(sample_table$karyotype))
          )
        
        #print the results out
        print(
          paste0(
            "Finished round ",
            as.character(i),
            ",",
            as.character(j),
            ". # of X responsive genes: ",
            as.character(numSig_npx),
            "; # of PAR-responsive genes: ",
            as.character(numSig_par)
          )
        )
        
        
        #advance j
        j <- j + 1
      }
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
colnames(saturationResults) <-
  c(
    "round",
    "npx_sig",
    "npx_sig_up",
    "npx_sig_up_ypair",
    "npx_sig_up_noY",
    "npx_sig_down",
    "par_sig",
    "par_sig_up",
    "par_sig_down",
    names(summary(sample_table$karyotype))
  )
write.table(
  saturationResults,
  file = paste0("saturationResults_100iterations_", as.character(i), ".txt"),
  quote = FALSE,
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t"
)

sessionInfo()
