library(DESeq2)
setwd("~/TCGA/")
list.files()

countData <- read.delim("Organized_Counts.csv", row.names = 1, header = TRUE, sep = ",")
head(countData)
groups <-read.csv("Metadata.csv", row.names = NULL, header = T, sep = ",")
head(groups)

#Filtering out genes w/ below 10 counts avg:
keep <- countData[rowMeans(countData) >= 10,]

# Making the dataset

dds <- DESeqDataSetFromMatrix(keep, colData = Metadata , design = ~ WFDC2)

dds <- DESeq(dds)
res <- results(dds) 

#Make the res into a data frame:
resdata<- as.data.frame(res)

#Convert row names into first column:
resdata <- cbind(rownames(resdata), data.frame(resdata, row.names=NULL))

#Assigning gene names reference file to human reference:
ref <- read.csv("~/genome_data/Human_Gene_Reference_GRCh38.csv", row.names=NULL)

#Merging resdata and ref:
merged <- as.data.frame(merge(resdata, ref, by.x = "rownames(resdata)",by.y = "Gene_ID", all.x = T))
merged2 <- merged[, c(8:16, 1:7)]
colnames(merged2) <- c("Chromosome", "Database", "ID_Type", "Start", "End","Strand","Gene_Version", "Gene_Name", "Gene_Type", "Gene_ID", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")

write.csv(merged2, row.names= FALSE, file = "results_WFDC2_high_vs_low.csv")
