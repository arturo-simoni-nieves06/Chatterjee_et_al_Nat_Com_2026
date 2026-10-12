# Load required libraries
library(ggplot2)
library(dplyr)

# Create the new data frame
data <- data.frame(
  Pathway = c(
    "MHC I mediated antigen presentation",
    "Microautophagy Signaling Pathway",
    "TGF-β Signaling",
    "Ferroptosis Signaling Pathway",
    "NRF2-mediated Stress Response",
    "PTEN Signaling",
    "Glutathione-mediated Detoxification",
    "PPAR Signaling"
  ),
  Zscore = c(
    -4.315,
    -3.317,
    -2.333,
    -1.508,
    -1.134,
    1.069,
    1.342,
    2.121
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
  scale_x_continuous(limits = c(-5, 3), breaks = seq(-5, 3, 1)) + # Adjusted x-axis limits for the new data
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
