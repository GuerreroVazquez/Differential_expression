
library( DESeq2 )
# library(ggplot2)
library(tidyverse)
library(tximport)
FC <- 1 # Fold-change cutoff
FDR <- 0.1 # FDR cutoff
alpha <- 0.1 # independent filtering, default



all_dfs <- c("GSE152558", "GSE157585", "GSE164471", "GSE164471_MAvO", "GSE164471_MAvY")

file = paste("/home/karen/Documents/phd/DifferentialExpression/gencode.v28.annotation.tx2gene.csv")
tx2gene <- read.csv(file, header = TRUE)
for (experiment in all_dfs)
{
  abundances_paths=paste0("/home/karen/Documents/phd/DifferentialExpression/",experiment,"/abundances.csv")
  df <- read.csv(abundances_paths, header = FALSE)
  files <- c()
  for( i in df[1]) 
    {  files <- append(files, i)
    }
  
  # files <- file.path(files=dir, type="kallisto", samples$run, "abundance.tsv")
  # names(files) <- paste0("sample", 1:6)
  
      files <- sort(files)
      txi<- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = FALSE)
      eval_samples <- read.csv(paste0("/home/karen/Documents/phd/DifferentialExpression/",experiment,"/metadata_filtered.csv"),                                file = header = TRUE)
      eval_samples <-eval_samples[order(eval_samples$Run),]
      sampleTable <- data.frame(condition = factor(eval_samples$category))
      dds = DESeqDataSetFromTximport(txi, colData=sampleTable, design=~condition )
      dds <- DESeq(dds, test="Wald")
      res <- results(dds)
      vsdata <- vst(dds, blind=FALSE)
      dists <- dist(t(assay(vsdata)))
      #plotPCA(vsdata) #using the DESEQ2 plotPCA fxn we can
      #plotMA(dds)

      no_na= na.omit(res)
      
      res <- DESeq2::results(dds, 
        lfcThreshold = log2(FC),
        altHypothesis = "greaterAbs",
        independentFiltering = FALSE,
        alpha = alpha
      )
      # Get up and down regulated
      res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC))
      res_up_regulates <- subset(res, padj < FDR & log2FoldChange > log2(FC))
      res_down_regulates <- subset(res, padj < FDR & log2FoldChange < -log2(FC))
    
      write.csv(res_up_regulates, paste0("UP_",experiment,"def_expressed.csv"))
      write.csv(res_down_regulates, paste0("DOWN_",experiment,"def_expressed.csv"))
}
