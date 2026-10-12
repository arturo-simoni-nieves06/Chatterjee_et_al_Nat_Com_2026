# Load required libraries
library(pheatmap)
library(RColorBrewer)
library(scales)  # For color manipulation

# Define the dataframe
data <- data.frame(
  symbol = c("Glrx3", "Gpx1", "Prdx4", "Sod3", "Gpx4", "Sod1", "Sod2", "Prdx6", "Prdx2", "Hmox2", 
             "Cat", "Prdx1", "Glrx", "Prdx5", "Prdx3", "Gsta4", "Gstm2", "Mgst1", "Gstm1", 
             "Gsta2", "Gstt1", "Gstm4", "Gstm5"),
  WT = c(159.7332883, 1005.82726, 363.5169297, 0.203310177, 317.4481953, 2230.847641, 
         624.4271509, 1264.176569, 3543.538286, 381.180566, 199.4705216, 84.76864754, 
         166.0551608, 2249.187622, 918.8823178, 2644.805708, 1392.181773, 5609.116507, 
         2956.150617, 29.01351642, 15.75323661, 23.65647225, 30.06555709),
  `p53 (-/-)` = c(505.5919159, 1814.138446, 558.6749409, 2.435919612, 330.6936164, 
                  2468.730482, 496.9080985, 1273.798827, 1968.238281, 314.2940855,
                  104.1224342, 128.4878472, 138.7704597, 1093.761364, 501.1676885,
                  851.6224106, 102.7431585, 578.4917899, 34.27072902, 
                  2.454570107, 0, 0, 0),
  mut175 = c(201.3484905, 1591.493169, 553.2221035, 3.065671046, 345.1722876, 2702.850684,
             561.6305191, 1368.808714, 3316.615935, 404.7215068, 139.2113106, 86.97553971,
             193.3813327, 1617.325428, 638.3615406, 2715.805187, 354.7039052, 2699.952369,
             317.4410508, 6.58334753, 2.886822664, 0, 4.306022596),
  mut273 = c(247.7622497, 1591.464053, 609.6737837, 3.118530077, 575.1703165, 2985.80012,
             574.7749285, 1392.192851, 3447.927009, 408.0675622, 146.4674881, 95.05729778,
             211.8554854, 1800.276408, 611.5811616, 3696.673155, 717.9308443, 3277.055,
             642.1853169, 10.01648431, 2.835544319, 0.283185494, 4.442793514)
)

# Define gene families and NRF2 regulation
gene_data <- data.frame(
  Gene = data$symbol,
  Family = c("GLRX", "GPX", "PRDX", "SOD", "GPX", "SOD", "SOD", "PRDX", "PRDX", "HMOX", 
             "CAT", "PRDX", "GLRX", "PRDX", "PRDX", "GST", "GST", "GST", "GST", 
             "GST", "GST", "GST", "GST"),
  NRF2_regulated = c(rep("Yes", 22), "Unknown")
)

# Order the data by Family
gene_data <- gene_data[order(gene_data$Family), ]
data <- data[match(gene_data$Gene, data$symbol), ]

# Prepare the matrix for heatmap
matrix_data <- as.matrix(data[, -1])
rownames(matrix_data) <- data$symbol

# Create annotation dataframe
annotation_df <- data.frame(
  Family = gene_data$Family,
  NRF2_regulated = gene_data$NRF2_regulated,
  row.names = gene_data$Gene
)

# Define colors for annotations (Modify these colors as needed)
ann_colors <- list(
  Family = c(SOD = "#E41A1C", GPX = "#377EB8", PRDX = "#4DAF4A", CAT = "#984EA3", 
             GLRX = "#FF7F00", HMOX = "#FFFF33", GST = "#A65628"),
  NRF2_regulated = c(Yes = "#66C2A5", Unknown = "#FC8D62")
)

# Define heatmap colors (Modify these as needed)
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Generate heatmap
pheatmap(matrix_data,
         annotation_row = annotation_df,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = FALSE,  # This prevents clustering of rows
         cluster_cols = TRUE,
         scale = "row",  # Scale rows (genes)
         color = heatmap_colors,
         main = "Gene Expression Heatmap with Family and NRF2 Regulation",
         fontsize = 8,
         cellwidth = 15,
         cellheight = 12)

# If you want to save the heatmap, uncomment the following lines:
# png("heatmap.png", width = 800, height = 600)
# pheatmap(...)  # Copy the pheatmap function call from above
# dev.off()
