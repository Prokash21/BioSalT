

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)  
library(gridExtra)

result <- read.csv("eADMSC_resLFC_data_for_heatmap_bar.csv")
################################################################################


genes <- result$Gene
conditions <- c("eMCF7", "eHeLa", "eMDAMB231", "eA549", "eH1975")

data_long <- result %>%
  select(Gene, ends_with("log2FoldChange"), ends_with("lfcSE"), ends_with("padj")) %>%
  pivot_longer(cols = -Gene, names_to = "Condition", values_to = "Value") %>%
  separate(Condition, into = c("Condition", "Type"), sep = "_") %>%
  pivot_wider(names_from = "Type", values_from = "Value") %>%
  rename(Log2FC = log2FoldChange, SE = lfcSE, Padj = padj)

# Define significance levels
data_long$Significance <- cut(data_long$Padj, 
                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                              labels = c("***", "**", "*", ""))

# Plot
ggplot(data_long, aes(x = Gene, y = Log2FC, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = Log2FC - SE, ymax = Log2FC + SE), 
                position = position_dodge(width = 0.7), width = 0.2) +
  geom_text(aes(label = Significance, y = Log2FC + SE + 0.2),
            position = position_dodge(width = 0.7), size = 5) +
  theme_minimal() +
  labs(title = "Relative Log2 Fold Change with Error Bars", 
       y = "Log2 Fold Change", x = "Gene",
       caption = "Significance: *** (p < 0.001), ** (0.001 ≤ p < 0.01), * (0.01 ≤ p < 0.05)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


