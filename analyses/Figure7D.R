# Load required libraries
library(tidyverse)
library(reshape2)
library(viridis)
library(svglite)  # For high-quality SVG output

# Create the dataframe
df <- data.frame(
  Cell_Type = c(rep("MΦ", 7), rep("CTLs", 14), rep("Other", 10)),
  Cytokines = c("CCL3", "CCL5", "CCL2", "IL4", "IL6", "IL13", "GCSF",
                "CCL2", "CCL3", "CCL5", "IL17", "Flt3", "Eselectin", "Pselectin",
                "IL7", "IL12", "IL15", "IL2", "IL4", "IL6",
                "CXCL1", "Gas6", "EGF", "MMP9", "HGF", "Endoglin", 
                "GMCSF", "IFNg", "TNFa", "CD14", "C5A"),
  mAvelumab_1 = c(5.73493899, 3.13483976, 0.63489852, 0.74788227, 0.06969057,
                  0.81996643, 1.088276, 0.63489852, 5.73493899, 3.13483976,
                  0.57513992, 2.31071165, 8.55210622, 1.14292683, 1.12120664,
                  2.19608933, 1.26170381, 1.91383016, 0.74788227, 0.06969057,
                  1.47302, 2.35145769, 1.02166447, 3.55263689, 0.9921876,
                  0.8627469, 0.46121291, 1.08902963, 0.97893026, 0.87510485,
                  0.74466664),
  mAvelumab_2 = c(5.77567346, 3.1071056, 0.63229643, 0.73998079, 0.06880926,
                  0.82476843, 1.08536064, 0.63229643, 5.77567346, 3.1071056,
                  0.56806796, 2.29027781, 8.58306815, 1.14891958, 1.1156063,
                  2.2017151, 1.25863809, 1.91335302, 0.73998079, 0.06880926,
                  1.44880854, 2.37485825, 1.01059171, 3.53229641, 0.99173641,
                  0.86539104, 0.45391631, 1.07297722, 0.98063873, 0.86937285,
                  0.73644374),
  AR37_1 = c(9.64249544, 12.8174937, 2.32747559, 1.19506704, 0.03482626,
             0.57500358, 1.02962198, 2.32747559, 9.64249544, 12.8174937,
             0.07873821, 7.83084134, 3.68846449, 1.38762781, 1.39611167,
             2.48031647, 1.12314793, 2.6190982, 1.19506704, 0.03482626,
             5.18945398, 1.35999047, 12.5939033, 3.69159597, 0.67971577,
             0.66117412, 0.66535489, 0.34697076, 2.18009101, 1.37416292,
             0.1093667),
  AR37_2 = c(9.70966573, 12.7711464, 2.28729515, 1.19687601, 0.03468828,
             0.56696801, 1.0228662, 2.28729515, 9.70966573, 12.7711464,
             0.07788744, 7.84600742, 3.71543168, 1.38369309, 1.39992812,
             2.47753236, 1.12562553, 2.60131896, 1.19687601, 0.03468828,
             5.20724, 1.37210362, 12.5939811, 3.66843564, 0.68521776,
             0.65972097, 0.66970024, 0.34970352, 2.19770741, 1.37519582,
             0.1096859),
  mAvelumab_AR37_1 = c(10.65007098, 9.90652833, 0.43142935, 0.56107565,
                       0.00565515, 0.41116385, 0.51030179, 0.43142935,
                       1.65007098, 9.90652833, 0.02202672, 3.04511195,
                       3.68149342, 1.43016925, 1.34486547, 1.24393226,
                       1.37441503, 5.08373857, 0.56107565, 0.00565515,
                       0.53868818, 0.67047681, 0.49999127, 0.74826269,
                       0.45925879, 0.39490444, 0.28288835, 0.25213603,
                       0.79600117, 0.33827747, 0.0457282),
  mAvelumab_AR37_2 = c(10.63261966, 9.7407475, 0.42742255, 0.56135948,
                       0.00568614, 0.40857883, 0.50173528, 0.42742255,
                       1.63261966, 9.7407475, 0.02219507, 2.99248818,
                       3.61040259, 1.40163849, 1.31879349, 1.23831542,
                       1.34818739, 5.09048207, 0.56135948, 0.00568614,
                       0.53231389, 0.65995458, 0.49196797, 0.74212159,
                       0.46214956, 0.3905912, 0.28466915, 0.25203763,
                       0.78012263, 0.33921573, 0.04495665)
)

# Calculate means for each treatment
df_mean <- df %>%
  mutate(
    mAvelumab = (mAvelumab_1 + mAvelumab_2)/2,
    AR37 = (AR37_1 + AR37_2)/2,
    mAvelumab_AR37 = (mAvelumab_AR37_1 + mAvelumab_AR37_2)/2
  ) %>%
  select(Cell_Type, Cytokines, mAvelumab, AR37, mAvelumab_AR37)

# Reshape data for heatmap
df_melt <- melt(df_mean, id.vars = c("Cell_Type", "Cytokines"),
                variable.name = "Treatment", value.name = "Value")

# Apply z-score normalization within each cytokine
df_melt_z <- df_melt %>%
  group_by(Cytokines) %>%
  mutate(
    Value_z = (Value - mean(Value)) / sd(Value),
    Value_scaled = (Value - min(Value)) / (max(Value) - min(Value))
  ) %>%
  ungroup()

# Create ordering factor to maintain original order
cytokine_order <- unique(df_melt_z$Cytokines)
df_melt_z$Cytokines <- factor(df_melt_z$Cytokines, levels = cytokine_order)

# NOW create the heatmap
heatmap_z <- ggplot(df_melt_z, aes(x = Treatment, y = Cytokines)) +
  geom_tile(aes(fill = Value_z), color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(Value, 2)), size = 3.2,
            color = ifelse(abs(df_melt_z$Value_z) > 1, "white", "black"),
            fontface = "bold") +
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "#ffffbf",
    high = "#d7191c",
    midpoint = 0,
    name = "Z-score",
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  facet_grid(Cell_Type ~ .,
             scales = "free_y",
             space = "free_y",
             switch = "y") +
  labs(
    title = "Cytokine Expression Heatmap (Z-score normalized)",
    subtitle = "Values in squares are original expression levels",
    x = "Treatment Condition",
    y = "Cytokines",
    caption = "Blue: Below average expression\nRed: Above average expression"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "gray40", hjust = 0),
    strip.text.y = element_text(angle = 0, size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = "gray70"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.3, "cm"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

# Display the heatmap
print(heatmap_z)

# SAVE AS SVG FOR INKSCAPE EDITING
ggsave("cytokine_heatmap_zscore.svg", heatmap_z,
       width = 10, height = 12, dpi = 800, bg = "white")

# Also save as high-resolution PNG for reference
ggsave("cytokine_heatmap_zscore.png", heatmap_z,
       width = 10, height = 12, dpi = 800, bg = "white")

cat("Heatmap saved as 'cytokine_heatmap_zscore.svg' - ready for editing in Inkscape!\n")
cat("\nTo edit in Inkscape:\n")
cat("1. Open the SVG file in Inkscape\n")
cat("2. Ungroup elements (Ctrl+Shift+G) to edit individually\n")
cat("3. Double-click text to modify\n")
cat("4. Select elements to change colors\n")