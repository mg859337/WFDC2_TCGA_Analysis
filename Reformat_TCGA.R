library(dplyr)

# First, figure out the median WFDC2 FPKM

FPKM <- read.csv("~/TCGA_FPKM.csv", row.names=NULL)
head(FPKM)
WFDC2 <- filter(FPKM, ensembl_gene_id == "ENSG00000101443")
WFDC2_Table <- as.data.frame(t(WFDC2))
colnames(WFDC2_Table) <- WFDC2_Table[1,]
write.csv(WFDC2_Table, "WFDC2_Table.csv", row.names = T, quote = F)

# Next, rearrange the count table so that we can easily make a metadata file that matches

WFDC2_Table <- read.csv("~/WFDC2_Table.csv", row.names=1)

median(WFDC2_Table$ENSG00000101443)
# 798.9293

high <- filter(WFDC2_Table, ENSG00000101443 > 798.9293)
low <- filter(WFDC2_Table, ENSG00000101443 < 798.9293)

highnames <- rownames(high)
lownames <- rownames(low)

write.csv(highnames, "High_WFDC2_Sample_Names.csv", row.names=F, quote=F)
write.csv(lownames, "Low_WFDC2_Sample_Names.csv", row.names=F, quote=F)

htseq_tcga_count_table <- read.csv("~/TCGA/htseq_tcga_count_table.csv", row.names=1)

high_counts <- select(htseq_tcga_count_table, highnames)
head(high_counts)
low_counts <- select(htseq_tcga_count_table, lownames)
head(low_counts)

organized <- merge(high_counts, low_counts, by= "row.names")
head(organized)

write.csv(organized, "Organized_Counts.csv", row.names = F, quote=F)
