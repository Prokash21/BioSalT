#install.packages("ape")
#install.packages("ggplot2")
#install.packages("ggtree")

#BiocManager::install("ape")
#BiocManager::install("ggtree")

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
getwd()

#PS
#set the working directory
setwd("E:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL/June_Run/GSE233233_QZ_006")


####################################### UMAP #######################################

# Load libraries
library(umap)
library(ggplot2)


# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data1 <- read.csv("Count_Data_final.csv", row.names = 1)
gene_expression_data2 <- gene_expression_data1 [1:18]
# Transpose the gene expression data
gene_expression_data <- t(gene_expression_data2)

# Example metadata (assuming you are reading from CSV files)
metadata1 <- read.csv("Meta_Data_final.csv")
metadata <- metadata1[1:18,]
library(ape)
library(ggplot2)
library(ggtree)

set.seed(123)
gene_expression_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(gene_expression_data) <- paste("Gene", 1:10, sep = "")
colnames(gene_expression_data) <- paste("Sample", 1:10, sep = "")
gene_expression_data <- gene_expression_data2
# Compute distance matrix
dist_matrix <- dist(t(gene_expression_data))

# Perform hierarchical clustering
hc <- hclust(dist_matrix)

# Convert to a phylogenetic tree
phylo_tree <- as.phylo(hc)


# Define groups for coloring (example)
groups <- rep(1:2, each = 5)
groups <- rep(metadata$agent)

# Create a data frame for group information
group_df <- data.frame(label = phylo_tree$tip.label, group = factor(groups))

# Visualize the phylogenetic tree with colors
ggtree(phylo_tree) %<+% group_df +
  geom_tiplab(aes(color = group), size = 3) +
  scale_color_manual(values = c("red", "blue")) +
  theme(legend.position = "right")

################################################## outliers ########################

# Load libraries
library(umap)
library(ggplot2)


# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data1 <- read.csv("Normalized_Salt_vs_control.csv", row.names = 1)
gene_expression_data1 <- read.csv("Remove_4_Normalized_Salt_vs_control.csv", row.names = 1)

gene_expression_data2 <- gene_expression_data1 [1:18]
gene_expression_data2 <- gene_expression_data2[, -c(2, 4, 5, 6)]

# Transpose the gene expression data
gene_expression_data <- t(gene_expression_data2)

# Example metadata (assuming you are reading from CSV files)
metadata1 <- read.csv("Meta_Data_final.csv")
metadata <- metadata1[1:18,]
metadata <- metadata [-c(2, 4, 5, 6),]
library(ape)
library(ggplot2)
library(ggtree)

set.seed(123)
gene_expression_data <- matrix(rnorm(100), nrow = 10, ncol = 10)
rownames(gene_expression_data) <- paste("Gene", 1:10, sep = "")
colnames(gene_expression_data) <- paste("Sample", 1:10, sep = "")
gene_expression_data <- gene_expression_data2
# Compute distance matrix
dist_matrix <- dist(t(gene_expression_data))

# Perform hierarchical clustering
hc <- hclust(dist_matrix)

# Convert to a phylogenetic tree
phylo_tree <- as.phylo(hc)


# Define groups for coloring (example)
groups <- rep(1:2, each = 5)
groups <- rep(metadata$agent)

# Create a data frame for group information
group_df <- data.frame(label = phylo_tree$tip.label, group = factor(groups))

# Visualize the phylogenetic tree with colors
ggtree(phylo_tree) %<+% group_df +
  geom_tiplab(aes(color = group), size = 3) +
  scale_color_manual(values = c("red", "blue")) +
  theme(legend.position = "right")

