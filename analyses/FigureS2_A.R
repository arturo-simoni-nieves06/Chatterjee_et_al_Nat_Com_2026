# ==========================================================
# 1. Load Libraries
# ==========================================================
library(ggplot2)
library(dplyr)
library(DESeq2) # Used for VST normalization (standard for bulk RNA-seq)

# ==========================================================
# 2. Load Data (Replace this with your real data)
# ==========================================================

# Create a count matrix (1000 genes, 12 samples)
counts <- matrix(rnbinom(n = 1000 * 12, mu = 100, size = 1/0.5), ncol = 12)
rownames(counts) <- paste0("Gene", 1:1000)
colnames(counts) <- paste0("Sample", 1:12)

# Create sample metadata matching the figure groups
# Based on the plot, we have 6 distinct groups
coldata <- data.frame(
  row.names = colnames(counts),
  Condition = factor(rep(c("WT", "p53(-/-) R270H", "p53(-/-) R172H", 
                           "PD-L1 KO", "p53(-/-) PD-L1 KO", "p53(-/-)"), each = 2))
)

# ==========================================================
# 3. Normalization (VST)
# ==========================================================
# PCA should be performed on transformed data, not raw counts
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ Condition)

# Run Variance Stabilizing Transformation
vsd <- vst(dds, blind = TRUE) # blind=TRUE is standard for QC/PCA

# ==========================================================
# 4. Calculate PCA
# ==========================================================
# Transpose the data so samples are rows and genes are columns
pca_data <- assay(vsd) %>% t() 

# Run PCA
pca_res <- prcomp(pca_data, scale. = FALSE)

# Extract coordinates for plotting
pca_df <- as.data.frame(pca_res$x)
pca_df$Condition <- coldata$Condition

# Calculate percentage variance for axis labels
percent_var <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)), 2)

# ==========================================================
# 5. Custom Styling and Plotting
# ==========================================================

# Define custom colors matching the image
custom_colors <- c(
  "WT" = "#D98880",               # Salmon/Red
  "p53(-/-) R270H" = "#52BE80",   # Green
  "p53(-/-) R172H" = "#BB8FCE",   # Purple
  "PD-L1 KO" = "#5D6D7E",         # Dark Blue/Grey
  "p53(-/-) PD-L1 KO" = "#C0392B",# Dark Red
  "p53(-/-)" = "#A2D9CE"          # Light Teal
)

# Create the plot
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
  # Points
  geom_point(size = 4, alpha = 0.9) +
  
  # Custom colors
  scale_color_manual(values = custom_colors) +
  
  # Axis Labels (Dynamic based on calculation)
  labs(
    title = "PCA plot: ID8 isogenic lines",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  
  # Theme customization to match the image style
  theme_bw() + # White background with border
  theme(
    # Title styling
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    
    # Axis text styling
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14, color = "black"),
    
    # Remove the legend (since the figure labels groups directly on the plot)
    legend.position = "none",
    
    # Grid styling
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  # Add text labels directly to the plot (Manual annotation)
  # Note: In a real analysis, you might use geom_text_repel from ggrepel 
  # to avoid overlapping, but here I am manually placing them to match the image.
  annotate("text", x = -22, y = -9, label = "PD-L1 KO", 
           fontface = "bold", size = 5, hjust = 0.5) +
  
  annotate("text", x = -15, y = 15, label = "p53 (-/-)\nR270H", 
           fontface = "bold", size = 5, hjust = 0.5) +
  
  annotate("text", x = -10, y = -3, label = "WT", 
           fontface = "bold", size = 5, hjust = 0.5) +
  
  annotate("text", x = 0, y = 17, label = "p53 (-/-)\nR172H", 
           fontface = "bold", size = 5, hjust = 0.5) +
  
  annotate("text", x = 15, y = 5, label = "p53 (-/-)\nPD-L1 KO", 
           fontface = "bold", size = 5, hjust = 0.5) +
  
  annotate("text", x = 25, y = -12, label = "p53 (-/-)", 
           fontface = "bold", size = 5, hjust = 0.5)

# Display plot
print(p)