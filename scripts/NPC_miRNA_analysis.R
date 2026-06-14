# =================================================
# Differential miRNA Expression Analysis in
# Nasopharyngeal Carcinoma (GSE70970)
# Author: Lopamudra Basu
# =================================================
# Required packages:
# install.packages("BiocManager")
# BiocManager::install("GEOquery")
# BiocManager::install("limma")
#=================================================
#Install GEOquery and limma
library(GEOquery)
library(limma)
# =================================================
#Download GEO Dataset
gset <- getGEO(
  "GSE70970",
  GSEMatrix = TRUE
)
#Extract Expression Matrix
expr <- exprs(gset[[1]])

#Inspect Metadata
pdata <- pData(gset[[1]])


expr_filtered <- expr[
  
  apply(expr, 1, var) != 0,
  
]
# =================================================

# =================================================
# Principal Component Analysis (PCA)

pca <- prcomp(
  t(expr_filtered),
  scale = TRUE
)

# Create sample groups
group <- factor(
  ifelse(
    pdata$source_name_ch1 ==
      "FFPE Nasopharyngeal Carcinoma Biopsy",
    "Tumor",
    "Normal"
  )
)

# Assign colors
colors <- ifelse(
  group == "Tumor",
  "red",
  "blue"
)

# Save PCA plot
png(
  "PCA_plot.png",
  width = 1200,
  height = 900
)

plot(
  pca$x[,1],
  pca$x[,2],
  col = colors,
  pch = 19,
  xlab = "PC1",
  ylab = "PC2",
  main = "Nasopharyngeal Carcinoma PCA"
)

legend(
  "topright",
  legend = c("Tumor", "Normal"),
  col = c("red", "blue"),
  pch = 19,
  cex = 0.8
)

dev.off()
# =================================================



# =================================================
#Design Matrix
design <- model.matrix(~ group)

#Differential Expression with limma
##Fit Linear Model
fit <- lmFit(expr_filtered, design)
###Empirical Bayes Moderation
fit <- eBayes(fit)
####Generate DEG Results
deg_results <- topTable(
  fit,
  
  coef = "groupTumor",
  
  number = Inf
)

#Significant miRNA
sig_deg <- deg_results[
  abs(deg_results$logFC) > 1 &
    deg_results$adj.P.Val < 0.05,
]
write.csv(
  deg_results,
  "Complete_DEG_Results.csv"
)

# =================================================


# =================================================
#create volcano plot

png("Volcano_plot.png", width=1200, height=900)

plot(
  deg_results$logFC,
  -log10(deg_results$adj.P.Val),
  
  pch = 16,
  
  col = ifelse(
    deg_results$adj.P.Val < 0.05 &
      abs(deg_results$logFC) > 1,
    "red",
    "grey"
  ),
  
  xlab = "log2 Fold Change",
  
  ylab = "-log10(FDR)",
  
  main = "Tumor vs Normal Volcano Plot"
)
dev.off()

# =================================================


# =================================================

#create heatmap of top 20 DEGs

png("Heatmap.png", width=1200, height=1000)
top20 <- rownames(
  sig_deg[order(sig_deg$adj.P.Val), ]
)[1:20]
heatmap_data <- expr_filtered[top20, ]
write.csv(
  heatmap_data,
  "Top20_Heatmap_miRNAs.csv"
)

heatmap(
  heatmap_data,
  
  scale = "row",
  
  col = heat.colors(100)
)
dev.off()

# =================================================


# =================================================

up_miRNA <- sig_deg[sig_deg$logFC > 1, ]
down_miRNA <- sig_deg[sig_deg$logFC < -1, ]


top_up <- head(up_miRNA,20)
top_down <- head(down_miRNA,20)

write.csv(top_up,
          "Top_Upregulated_miRNA.csv")

write.csv(top_down,
          "Top_Downregulated_miRNA.csv")

#Create summary table 
summary_table <- rbind(
  head(up_miRNA,10),
  head(down_miRNA,10)
)


write.csv(
  summary_table,
  "Top20_miRNA_Summary.csv"
)
#Top up-regulated and down-regulated miRNA 
png(
  "Top_miRNA_Barplots.png",
  width=1600,
  height=800,
  res=150
)
par(mfrow=c(1,2))

barplot(
  top_up$logFC,
  names.arg=rownames(top_up),
  las=2,
  col="red",
  main="Top Upregulated miRNAs"
)

barplot(
  abs(top_down$logFC),
  names.arg=rownames(top_down),
  las=2,
  col="blue",
  main="Top Downregulated miRNAs"
)

par(mfrow=c(1,1))

dev.off()

write.csv(
  sig_deg,
  "All_Significant_miRNAs.csv"
)

human_miRNA <- sig_deg[grep("^hsa", rownames(sig_deg)), ]

write.csv(
  human_miRNA,
  "Human_miRNA_Only.csv"
)

ebv_miRNA <- sig_deg[grep("^ebv", rownames(sig_deg)), ]

write.csv(
  ebv_miRNA,
  "EBV_miRNA_Only.csv"
)


cat(
  "Total significant miRNAs:",
  nrow(sig_deg),
  "\n"
)

cat(
  "Upregulated miRNAs:",
  nrow(up_miRNA),
  "\n"
)

cat(
  "Downregulated miRNAs:",
  nrow(down_miRNA),
  "\n"
)