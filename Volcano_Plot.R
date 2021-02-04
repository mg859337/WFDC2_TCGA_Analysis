library(ggplot2)
library(dplyr)
library(ggrepel)

# First load the data

results <- read.csv("~/results_WFDC2_high_vs_low.csv")

# Second, make sure we're specifying the DEGs 

cleaned <- filter(results, results$padj != "NA")
cleaned2 <- filter(cleaned, cleaned$Gene_Type == "protein_coding")
cleaned2$Significant <- ifelse(cleaned2$padj < 0.05 & cleaned2$log2FoldChange >= 0.5 | cleaned2$padj < 0.05 & cleaned2$log2FoldChange <= -0.5, "p-adj. < 0.05", "Not Sig")


# Third, make the plot

plot <- ggplot(cleaned2, aes(x=log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("gray50", "red")) +
  theme(legend.position = "bottom") 

plot

# Fourth, add the labels

target <- c("WFDC2","SLPI","LAMTOR4","CST3","SPINK1")

new <- filter(cleaned, Gene_Name %in% target)

plot + geom_text_repel(
  data = new,
  aes(label = Gene_Name),
  size = 5,
  box.padding = unit(0.45, "lines"),
  point.padding = unit(0.4, "lines"),
  fontface = 2
)
