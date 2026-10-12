# Load required libraries
library(ggplot2)
library(dplyr)

# Create the data frame
data <- data.frame(
  Pathway = c("ERK/MAPK Signaling", "Signaling by EGFR", "Actin Cytoskeleton Signaling", 
              "TGF-β Signaling", "Regulation of the EMT by Growth Factors", 
              "NO and ROS production in MΦ", "Netrin Signaling", 
              "Cell Junction Organisation", "Glutathione-mediated Detoxification", 
              "PPAR Signaling", "Fatty Acid α-oxidation", "MHC class II antigen presentation"),
  Zscore = c(2.673, 2, 1.897, 1.633, 1.606, 1.342, -1.069, -1.134, -1.414, -1.633, -2, -2.333)
)

# Order the data by Z-score
data <- data %>% arrange(Zscore)

# Create the lollipop plot
ggplot(data, aes(x = Zscore, y = reorder(Pathway, Zscore))) +
  geom_segment(aes(x = 0, xend = Zscore, y = Pathway, yend = Pathway), 
               color = "grey50") +
  geom_point(aes(color = Zscore > 0), size = 6) +
  scale_color_manual(values = c("blue", "red")) +
  labs(x = "Z-score", y = NULL, title = "Pathway Enrichment Analysis") +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_text(size = 10)
  ) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
