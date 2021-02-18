#!/usr/bin/env Rscript

# Rscript for BCR venn diagram ploting
# (based on Wenwen Xiang @ TakaraBio script)
# TakaraBio USA
# Author: Wenwen Xiang
# version0.1
# Data: May 4th 2018
# Usage: Rscript bcr_venn.R folder_dirctory report_name
# use  median+3*sd as negctrl filter

library(dplyr)
library(limma)

# --------------------------------------------
# input: .aa_count.txt file with examples like
# cloneCount,aaSeqCDR3
# 96.0,CARSPIDSSYLFDYW
# 88.0,CAKDRSVGTVVVGRYFDYW
# 76.0,CARDDDLEYVRADLDFW
# 73.0,CARVQKNYYYGMDVW
# --------------------------------------------

# --------------------------------------------
#     define variables & set analysis directory
# --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
report_name <- args[1]
i <- ((length(args)+1)/2)
files <- args[2:i]
samples <- args[(i+1):length(args)]
venn_colors <- c("lightblue", "purple", "orange", "lightgreen", "grey")

if (length(files)<2) {
  print("Need 2 or more files to draw Venn Diagram. Quit.");
  quit(status=0);
}
if (length(files)!=length(samples)) {
  print("Number of files provided don't match number of sample names provided. Quit.");
  quit(status=0);
}

counts <- list();
unique_cdr3_seqs <- c();


for (file in files) {
  x <- read.table(file, sep=',', header=T, stringsAsFactors=F);
  if (length(grep("cloneCount", colnames(x)))==0) {
    print(paste("ERROR: Cannot find column 'cloneCount' in ", file, "...abort!",sep=''));
    quit(status=1);
  }
  if (length(grep("aaSeqCDR3", colnames(x)))==0) {
    print(paste("ERROR: Cannot find column 'aaSeqCDR3' in ", file, "...abort!",sep=''));
    quit(status=1);
  }

  x.summed <- x %>% group_by(aaSeqCDR3) %>% summarise(aa_count=sum(cloneCount))
  #colnames(x)[grep("cloneCount", colnames(x))] <- paste("cloneCount", file, sep='.')
  counts <- append(counts, list(x.summed))
  unique_cdr3_seqs <- union(unique_cdr3_seqs, x.summed$aaSeqCDR3);
}

df <- data.frame(aaSeqCDR3=unique_cdr3_seqs)
for (i in 1:length(counts)) {
  df[match(counts[[i]]$aaSeqCDR3, unique_cdr3_seqs),paste("aa_count",samples[i],sep='.')] <- counts[[i]]$aa_count;
}

write.table(df,file=paste(report_name,"_ordered_clono_count.csv", sep=""),sep=",",quote=F,col.names=T,row.names=F)

m <- as.matrix(df[,2:dim(df)[2]])
m[is.na(m)] <- 0
venn_counts <- vennCounts(m)
png(paste(report_name,".png",sep=""), width=3.25, height=3.25, units="in", res=1200, pointsize=8)
vennDiagram(venn_counts, main=paste(report_name," Venn Diagram", sep=""), names=samples, circle.col=venn_colors[1:length(samples)])
dev.off()

