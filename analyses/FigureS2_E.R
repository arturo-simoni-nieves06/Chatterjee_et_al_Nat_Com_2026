# Load required libraries
library(ggplot2)
library(dplyr)

# Create the data frame with your new data
data <- data.frame(
  Pathway = c(
    "Macropinocytosis Signaling",
    "Autophagy",
    "Ferroptosis Signaling Pathway",
    "eNOS Signaling",
    "Mitochondrial Dysfunction",
    "Glutathione-mediated Detoxification",
    "Fatty Acid α-oxidation",
    "Interferon Signaling"
  ),
  Zscore = c(
    -2.236,
    -1.528,
    -1.387,
    -1.134,
    1.279,
    1.89,
    2,
    2.449
  )
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
