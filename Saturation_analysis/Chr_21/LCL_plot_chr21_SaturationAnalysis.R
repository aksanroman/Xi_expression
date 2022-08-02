#Making plots

suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("ggplot2"))

set.seed(seed = 1)

#read in all of the results
results45 <- data.frame(Sample_num = rep(45),read.delim(file="saturationResults_100iterations_45.txt"))
results10 <- data.frame(Sample_num = rep(10),read.delim(file="saturationResults_100iterations_10.txt"))
colnames(results10) <- colnames(results45)
results15 <- data.frame(Sample_num = rep(15),read.delim(file="saturationResults_100iterations_25.txt"))
colnames(results15) <- colnames(results45)
results20 <- data.frame(Sample_num = rep(20),read.delim(file="saturationResults_100iterations_20.txt"))
colnames(results20) <- colnames(results45)
results25 <- data.frame(Sample_num = rep(25),read.delim(file="saturationResults_100iterations_25.txt"))
colnames(results25) <- colnames(results45)
results30 <- data.frame(Sample_num = rep(30),read.delim(file="saturationResults_100iterations_30.txt"))
colnames(results30) <- colnames(results45)
results35 <- data.frame(Sample_num = rep(35),read.delim(file="saturationResults_100iterations_35.txt"))
colnames(results35) <- colnames(results45)
results40 <- data.frame(Sample_num = rep(40),read.delim(file="saturationResults_100iterations_40.txt"))
colnames(results40) <- colnames(results45)


#make the table
saturationResults <- NULL
saturationResults <- rbind(results10, results15, results20, results25, results30, results35, results40, results45)

#Plot the results to see whether it is saturating. 
data_summary <- function(x){
  # x <- gnomad_pLI_curated_NPX
  m <- median(x, na.rm=TRUE)
  ymin <- quantile(x, na.rm=TRUE)[2]
  ymax <- quantile(x, na.rm=TRUE)[4]
  test <- c("y"=m,"ymin"=ymin,"ymax"=ymax)
  names(test) <- c("y","ymin","ymax")
  return(test)
}

pdf(file = "saturation_plots_chr21_100iterations.pdf", width=2.5, height=2.5)

ggplot(data=saturationResults, mapping = aes(x=as.factor(Sample_num), y=chr21)) + 
  geom_violin(size=0.25, scale = "width", col="#00000000", fill="00000040") +
  stat_summary(fun.data=data_summary, size=0.25) +
  theme_classic(base_size = 8) + 
  theme(
    axis.text = element_text(color = "black"), 
    axis.ticks = element_line(color = "black", size = 0.25), 
    axis.line = element_line(size = 0.25),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(hjust = 0.5),
    legend.position = "none",
  ) +
  labs(
    title="Significant Chr21 genes, LCLs",
    y="Number of genes (FDR<0.05)",
    x="Sample size"
  ) 
dev.off()
