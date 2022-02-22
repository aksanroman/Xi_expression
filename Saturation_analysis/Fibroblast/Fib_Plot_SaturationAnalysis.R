#Making plots
suppressPackageStartupMessages(library("vioplot"))
suppressPackageStartupMessages(library("ggplot2"))

set.seed(seed = 1)
#read in all of the results
results10 <- data.frame(Sample_num = rep(10),read.delim(file="saturationResults_100iterations_10.txt"))
results20 <- data.frame(Sample_num = rep(20),read.delim(file="saturationResults_100iterations_20.txt"))
results30 <- data.frame(Sample_num = rep(30),read.delim(file="saturationResults_100iterations_30.txt"))
results40 <- data.frame(Sample_num = rep(40),read.delim(file="saturationResults_100iterations_40.txt"))
results50 <- data.frame(Sample_num = rep(50),read.delim(file="saturationResults_100iterations_50.txt"))
results60 <- data.frame(Sample_num = rep(60),read.delim(file="saturationResults_100iterations_60.txt"))
results70 <- data.frame(Sample_num = rep(70),read.delim(file="saturationResults_100iterations_70.txt"))
results80 <- data.frame(Sample_num = rep(80),read.delim(file="saturationResults_100iterations_80.txt"))
results90 <- data.frame(Sample_num = rep(90),read.delim(file="saturationResults_100iterations_90.txt"))
results99 <- data.frame(Sample_num = rep(99),read.delim(file="saturationResults_100iterations_99.txt"))

results10_male <- data.frame(Sample_num = rep(10),read.delim(file="saturationResults_100iterations_male_10.txt"))
results15_male <- data.frame(Sample_num = rep(15),read.delim(file="saturationResults_100iterations_male_15.txt"))
results20_male <- data.frame(Sample_num = rep(20),read.delim(file="saturationResults_100iterations_male_20.txt"))
results25_male <- data.frame(Sample_num = rep(25),read.delim(file="saturationResults_100iterations_male_25.txt"))
results30_male <- data.frame(Sample_num = rep(30),read.delim(file="saturationResults_100iterations_male_30.txt"))
results35_male <- data.frame(Sample_num = rep(35),read.delim(file="saturationResults_100iterations_male_35.txt"))
results40_male <- data.frame(Sample_num = rep(40),read.delim(file="saturationResults_100iterations_male_40.txt"))
results45_male <- data.frame(Sample_num = rep(45),read.delim(file="saturationResults_100iterations_male_45.txt"))
results50_male <- data.frame(Sample_num = rep(50),read.delim(file="saturationResults_100iterations_male_50.txt"))
results52_male <- data.frame(Sample_num = rep(52),read.delim(file="saturationResults_100iterations_male_52.txt"))


#make the table
saturationResults <- NULL
saturationResults_male <- NULL
saturationResults <- rbind(saturationResults,results10,results20, results30,results40,results50, results60, results70, results80, results90, results99)
saturationResults_male <- rbind(saturationResults_male,results10_male,results15_male, results20_male,results25_male, results30_male,results35_male, results40_male,
                                results45_male,results50_male,results52_male)

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

pdf(file = "saturation_plots_NPX_PAR_NPY_100iterations_fib.pdf", width=2.5, height=2.5)
ggplot(data=saturationResults, mapping = aes(x=as.factor(Sample_num), y=npx_sig_up_ypair)) + 
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
                title="NPX-NPY genes with deltaEx > 0, Fibs",
                y="Number of genes (FDR<0.05)",
                x="Sample size"
        ) + ylim(0,10)

ggplot(data=saturationResults, mapping = aes(x=as.factor(Sample_num), y=npx_sig_up_noY)) + 
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
                title="NPX no Y genes with deltaEx > 0, Fibs",
                y="Number of genes (FDR<0.05)",
                x="Sample size"
        ) + ylim(0,60)

ggplot(data=saturationResults, mapping = aes(x=as.factor(Sample_num), y=npx_sig_down)) + 
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
                title="NPX genes with deltaEx < 0 , Fibs",
                y="Number of genes (FDR<0.05)",
                x="Sample size"
        ) + ylim(0,60)

ggplot(data=saturationResults, mapping = aes(x=as.factor(Sample_num), y=par_sig_up)) + 
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
                title="Significant PAR genes, Fibs",
                y="Number of genes (FDR<0.05)",
                x="Sample size"
        ) + ylim(0,10)

ggplot(data=saturationResults_male, mapping = aes(x=as.factor(Sample_num), y=npy_sig)) + 
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
                title="NPY genes with deltaEy > 0, Fibs",
                y="Number of genes (FDR<0.05)",
                x="Sample size"
        ) + ylim(0,11)
dev.off()

