# First convert symbols to entrez IDs

library(org.Hs.eg.db)

DEGs <- read.table("DEGs.txt") # This is a list of DEGs from DESeq2 saved as a text file

degs <- as.character(DEGs$V1)

conversion <- mapIds(org.Hs.eg.db, degs, 'ENTREZID', 'ENSEMBL')


# Load up the packages needed

library(DOSE)
library(clusterProfiler)
library(ggplot2)

# For over representation analysis, we need is a gene ID vector 

genes <- conversion[!is.na(conversion)]

ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "fdr",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

head(ego) 

# Try to remove some redundancy

bp2 <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun = min)
head(bp2) 

write.csv(bp2, "EnrichGO2_DEGs.csv")

# Make dotplot of 15 categories

p1 <- dotplot(bp2, showCategory=15, orderBy = "GeneRatio", font.size = 16)
p1
