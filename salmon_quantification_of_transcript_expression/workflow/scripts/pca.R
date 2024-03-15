## Principal components analysis using log2TPM values ##

suppressPackageStartupMessages(library("edgeR"))
library("magrittr")
library("stringr")
library(ggfortify)
library(ggrepel)
library(ggthemes)
library("stringi")
library(dplyr)
library(tidyverse)

## Get input file from standard input ------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) > 0, file.exists(args))
f_counts <- args

## open file from stdin ##
tpm <- read.table(args[1], header=TRUE)

## make first column as rownames ##
rownames(tpm) <- tpm[,1]
tpm <- tpm[,-1]

## transform colnames to same syntax ##
colnames(tpm) <- str_replace(colnames(tpm), "results.|.s.bam|_", "")

## convert to matrix, to allow for log2 transformation ##
tpm <- as.matrix(tpm)

## function to convert tpm to log2tpm ##
logTPM <- function(tpm, dividebyten=TRUE) {
  if(dividebyten) {
    logtpm <- log(tpm/10+1, 2)}
  else if(!dividebyten) {
    logtpm <- log(tpm+1, 2)}
  return(logtpm)
}

## convert tpm to log2tpms ##
logtpms<-logTPM(tpm, dividebyten = FALSE)

## transpose ##
xt = t(logtpms)

## generate sample names for grouping replicates ##
xt <- as.data.frame(xt)
groups <- str_replace(rownames(xt), "\\d$", "")

## add column with sample name to group replicates ##
xtl <- xt %>% add_column(Sample = groups)

## run PCA analysis ##
pca_res = prcomp(xt, center=T, scale.=F)

svg("PCA.svg")

## draw auto-PCA plot with color mappings for groups ##
print ( autoplot(pca_res, data=xtl, colour='Sample', legend.size=5) +  labs(title = "Principal Component Analysis using log2TPM values" ) + theme_bw() + geom_text_repel(aes(label=rownames(xt),color=Sample), show.legend = FALSE  ) +
        theme(plot.title=element_text(face="bold",hjust=0.5), legend.title = element_text(size=12), legend.text = element_text(size=12), axis.title=element_text(size=12))   )

dev.off()
