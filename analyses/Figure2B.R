# Load required libraries
library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(here)

# Load normalized counts data and ensure unique row names
normalized_counts <- read.csv(here::here("data", "raw-data", "WT_vs_mut273_DEG.csv"), row.names = NULL)
rownames(normalized_counts) <- make.unique(as.character(normalized_counts[, 1]))
normalized_counts <- normalized_counts[, -1]

# Inspect the first few rows of the data
head(normalized_counts)

# Round the counts to the nearest integer (if necessary)
normalized_counts <- round(normalized_counts)

# Check for negative values
if (any(normalized_counts < 0)) {
  stop("Count data contains negative values.")
}

# Define conditions for samples based on your dataset
conditions <- c(rep("WT", 3), rep("mut273", 3))

# Create colData with matching sample names and correct length
colData <- data.frame(condition = factor(conditions, levels = c("WT", "mut273")), 
                      row.names = colnames(normalized_counts))

# Check dimensions to ensure they match
if (ncol(normalized_counts) != nrow(colData)) {
  stop("Mismatch between colData and normalized_counts: ",
       nrow(colData), " rows in colData vs ", ncol(normalized_counts), " columns in normalized_counts.")
}

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = as.matrix(normalized_counts), 
                              colData = colData, 
                              design = ~ condition)

# Run DESeq analysis
dds <- DESeq(dds)

# Get results from the DESeq analysis (mut273 vs WT)
results <- results(dds, contrast = c("condition", "mut273", "WT"))

# Create a data frame for plotting
volcano_data <- as.data.frame(results)
volcano_data$gene <- rownames(volcano_data)
volcano_data <- na.omit(volcano_data)  # Remove NA values

# Rename columns for easier access
colnames(volcano_data)[colnames(volcano_data) == "log2FoldChange"] <- "log2FC"
colnames(volcano_data)[colnames(volcano_data) == "padj"] <- "qvalue"

# Define the specific genes to highlight
emt_genes <- c("Twist2", "Col6a2", "Mmp2", "Timp3", "Slit3", "Ncam1")
ros_genes <- c("Gstm4", "Gstm5", "Chac1", "Gsta3", "Gstk1", "Gstm1")

# Add a column for gene type (ROS, EMT, or Other)
volcano_data$GeneType <- ifelse(volcano_data$gene %in% ros_genes, "ROS",
                                ifelse(volcano_data$gene %in% emt_genes, "EMT", "Other"))

# Add a column for regulation status
volcano_data$RegulationStatus <- ifelse(volcano_data$log2FC > 0.5 & volcano_data$qvalue < 0.05, "Upregulated in mut273",
                                        ifelse(volcano_data$log2FC < -0.5 & volcano_data$qvalue < 0.05, "Downregulated in mut273", "Not Significant"))

# Combine EMT and ROS genes
highlight_genes <- c(emt_genes, ros_genes)

# Create the volcano plot
volcano_plot <- ggplot(volcano_data, aes(x = log2FC, y = -log10(qvalue))) +
  geom_point(aes(color = RegulationStatus), alpha = 0.5) +
  geom_point(data = subset(volcano_data, gene %in% highlight_genes), 
             aes(color = GeneType), size = 3) +
  scale_color_manual(values = c("Upregulated in mut273" = "darkorange", 
                                "Downregulated in mut273" = "blue", 
                                "Not Significant" = "grey",
                                "EMT" = "purple", 
                                "ROS" = "darkgreen")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change (mut273 vs WT)",
       y = "-Log10 Adjusted P-value",
       title = "Volcano Plot: mut273 vs WT DEGs",
       color = "Gene Type") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black") +
  
  # Label the specified EMT and ROS genes
  geom_label_repel(data = subset(volcano_data, gene %in% highlight_genes), 
                   aes(label = gene, color = GeneType), 
                   size = 3, box.padding = 0.75, max.overlaps = Inf, 
                   fill = "white", alpha = 0.8, label.size = 0.75)

# Save the volcano plot as an SVG file
ggsave(filename = here::here("figures", "volcano_plot_mut273_vs_WT.svg"),
       plot = volcano_plot, width = 10, height = 8)

# Optionally display the plot in R
print(volcano_plot)