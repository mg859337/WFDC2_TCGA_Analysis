library(DESeq2)
library(ggplot2)

setwd("~/TCGA/")
list.files()
countData <- read.csv("Organized_Counts.csv", row.names = 1, header = TRUE, sep = ",")
head(countData)
groups <-read.csv("Metadata.csv",  header = T, sep = ",")
head(groups)

#Filtering out genes w/ below 10 counts avg:
keep <- countData[rowMeans(countData) >= 10,]

# Formatting the data into a dds
dds <- DESeqDataSetFromMatrix(keep, colData = Metadata, design = ~ WFDC2)

# This is the transformation that the PCA plot will be using
vsd <- vst(dds, blind=TRUE)

# Making the PCA plot
pcaData <- plotPCA(vsd, intgroup="WFDC2", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=WFDC2)) +
  geom_point(size=4) +
  theme(legend.text=element_text(size=16)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()