#
#

library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
getwd()

#set the working directory
setwd("E:/DWCT/M_Analysis_Ajit_Sir/GSE152620/June_Run")

#load the count data
count_data <- read.csv("Count_Data_final.csv", header=TRUE,row.names = 1)
colnames(count_data)
head(count_data)

#load the sample info
sample_info <- read.csv("Meta_Data_final.csv", header = TRUE,row.names = 1)
colnames(sample_info)
head(sample_info)

#set factor levels
sample_info$treatment <- factor(sample_info$treatment)

# Convert non-integer values to integers in count data
count_data <- round(count_data)
head(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)
head(new_count_data)

unique(sample_info$treatment)
#unique(sample_info$Time_point)
#unique(sample_info$tissue)

dim(count_data)
dim(new_count_data)
dim(sample_info)

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~ treatment)

# Perform DESeq2 analysis
dds <- DESeq(dds)
head(dds)

#set the factor level
dds$treatment <- factor(dds$treatment, levels = c ("Control","Salinity","Heat","Salinity + heat")) 

#filter the genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#set the referene for the treatment factor
dds$treatment <- relevel(dds$treatment , ref = "Control")
dds$treatment

#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)

#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"Normalized_Salt_vs_control.csv")

#Identify available coefficient names
coeff_names <- resultsNames(dds)

#Print the coefficient names
print(coeff_names)

#[1] "Intercept"                    "agent_Salt.stress_vs_control"

resLFC <- lfcShrink(dds, coef ="treatment_Salinity_vs_Control"  , type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)
##################18th oct
resLFC_p_cut <- as.data.frame(resLFC[ resLFC$padj<0.05,])
write.csv(resLFC_p_cut, file = "resLFC_p_cut.csv")
#################

# Create a Volcano Plot
ggplot(resLFC,aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Scatter plot points with color-coded regulation
  geom_point(aes(color = ifelse(log2FoldChange > 1.0 & -log10(padj) > 1.3, "Upregulated",
                                ifelse(log2FoldChange < -1.0 & -log10(padj) > 1.3, "Downregulated", "Not Significant"))),
             size = 2.5, alpha = 0.5) +
  
  # Add horizontal dashed line
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  
  # Add vertical dashed lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  # Customize plot labels and add the header
  labs(
    title = "HeatMap of agent_Salt.stress_vs_control", # Add the header here
    x = "Log2 Fold Change",
    y = "-log10(padj)",
    color = "Regulation"
  ) +
  
  # Customize color palette for regulation categories
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  
  # Use a minimal theme for the plot
  theme_minimal()

########################################## padj #############################################

Upregulated <- resLFC[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05, ]
top_UP_hits1 <- Upregulated[order(Upregulated$padj ),]
write.csv(top_UP_hits1, file = "UP_Salt_stress_vs_control.csv", row.names = TRUE )

Downregulated <- resLFC[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05, ]
top_hits1 <- Downregulated[order(Downregulated$padj),]
write.csv(top_hits1, file = "DOWN_agent_Salt_vs_control.csv", row.names = TRUE)

######################################### UMAP #######################################

# Load libraries
library(umap)
library(ggplot2)

# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data <- read.csv("Normalized_Salt_vs_control.csv", row.names = 1)
gene_expression_data <- gene_expression_data [ 1:6]
# Transpose the gene expression data
gene_expression_data <- t(gene_expression_data)

# Example metadata (assuming you are reading from CSV files)
metadata <- read.csv("Meta_Data_final.csv")
metadata <- metadata [1:6,]
# Check dimensions
print(dim(gene_expression_data))  # Check dimensions of gene expression data
print(dim(metadata))              # Check dimensions of metadata

# Ensure the number of samples match
if (nrow(gene_expression_data) != nrow(metadata)) {
  stop("Number of samples in gene_expression_data and metadata do not match!")
}

# Convert Time to factor
metadata$treatment <- factor(metadata$treatment)

# Set random seed for reproducibility
set.seed(123)

# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Extract UMAP coordinates and combine with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = treatment)) +
  geom_point(size = 3) +
  labs(title = "GSE152620_UMAP",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()


# Plot using ggplot2 with stat_ellipse and labels
ggplot(umap_df, aes(x = X1, y = X2, color = treatment)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.2) +  # Add confidence ellipses
  geom_text(aes(label = rownames(umap_df)), size = 2, nudge_y = 0.1) +  # Add sample labels
  labs(title = "GSE152620_UMAP",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()


################################################## PCA ############################3

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
#library(CorLevelPlot)
library(gridExtra)

data <- gene_expression_data
#data1 <- read.csv("Count_Data_final.csv", header = TRUE, row.names = 1)
data <- as.data.frame(data)
#metadata <- read.csv("Meta_Data_final.csv", header = TRUE, row.names = 1) 
metadata <- metadata
# Quality Control: detect outlier genes
gsg <- goodSamplesGenes(data)
summary(gsg)

# Remove outlier genes
data <- data[gsg$goodGenes == TRUE,]

# Detect outlier samples using hierarchical clustering
htree <- hclust(dist(data), method = "average")
######
#######
# Plot hierarchical clustering with colors based on metadata groups
plot(htree, labels = metadata$Sample.Name, main = "Hierarchical Clustering Dendrogram", col=metadata$treatment)
rect.hclust(htree, k = length(unique(metadata$group)), border = 2:5)


####################################### 29 ######################################


#set the working directory
setwd("E:/DWCT/M_Analysis_Ajit_Sir/GSE152620/June_Run")

#load the count data
count_data <- read.csv("Count_Data_final.csv", header=TRUE,row.names = 1)
count_data <- count_data[-5]
colnames(count_data)
head(count_data)

#load the sample info
sample_info <- read.csv("Meta_Data_final.csv", header = TRUE,row.names = 1)
sample_info <- sample_info[-5,]
colnames(sample_info)
head(sample_info)

#set factor levels
sample_info$treatment <- factor(sample_info$treatment)

# Convert non-integer values to integers in count data
count_data <- round(count_data)
head(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)
head(new_count_data)

unique(sample_info$treatment)
#unique(sample_info$Time_point)
#unique(sample_info$tissue)

dim(count_data)
dim(new_count_data)
dim(sample_info)

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~ treatment)

# Perform DESeq2 analysis
dds <- DESeq(dds)
head(dds)

#set the factor level
dds$treatment <- factor(dds$treatment, levels = c ("Control","Salinity","Heat","Salinity + heat")) 

#filter the genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#set the referene for the treatment factor
dds$treatment <- relevel(dds$treatment , ref = "Control")
dds$treatment

#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)

#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"29_X_Normalized_Salt_vs_control.csv")

#Identify available coefficient names
coeff_names <- resultsNames(dds)

#Print the coefficient names
print(coeff_names)

#[1] "Intercept"                    "agent_Salt.stress_vs_control"

resLFC <- lfcShrink(dds, coef ="treatment_Salinity_vs_Control"  , type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)
##################18th oct
resLFC_p_cut <- as.data.frame(resLFC[ resLFC$padj<0.05,])
write.csv(resLFC_p_cut, file = "29_X_resLFC_p_cut.csv")
#################

# Create a Volcano Plot
ggplot(resLFC,aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Scatter plot points with color-coded regulation
  geom_point(aes(color = ifelse(log2FoldChange > 1.0 & -log10(padj) > 1.3, "Upregulated",
                                ifelse(log2FoldChange < -1.0 & -log10(padj) > 1.3, "Downregulated", "Not Significant"))),
             size = 2.5, alpha = 0.5) +
  
  # Add horizontal dashed line
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  
  # Add vertical dashed lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  
  # Customize plot labels and add the header
  labs(
    title = "HeatMap of agent_Salt.stress_vs_control", # Add the header here
    x = "Log2 Fold Change",
    y = "-log10(padj)",
    color = "Regulation"
  ) +
  
  # Customize color palette for regulation categories
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  
  # Use a minimal theme for the plot
  theme_minimal()

########################################## padj #############################################

Upregulated <- resLFC[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05, ]
top_UP_hits1 <- Upregulated[order(Upregulated$padj ),]
write.csv(top_UP_hits1, file = "29_X_UP_Salt_stress_vs_control.csv", row.names = TRUE )

Downregulated <- resLFC[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05, ]
top_hits1 <- Downregulated[order(Downregulated$padj),]
write.csv(top_hits1, file = "29_X_DOWN_agent_Salt_vs_control.csv", row.names = TRUE)

############################################### umap #################################

# Load libraries
library(umap)
library(ggplot2)

# Example gene expression data (assuming you are reading from CSV files)
gene_expression_data <- count_data [ 1:5]
# Transpose the gene expression data
gene_expression_data <- t(gene_expression_data)

# Example metadata (assuming you are reading from CSV files)
#metadata <- read.csv("Meta_Data_final.csv")
metadata <- sample_info [1:5,]
# Check dimensions
print(dim(gene_expression_data))  # Check dimensions of gene expression data
print(dim(metadata))              # Check dimensions of metadata

# Ensure the number of samples match
if (nrow(gene_expression_data) != nrow(metadata)) {
  stop("Number of samples in gene_expression_data and metadata do not match!")
}

# Convert Time to factor
metadata$treatment <- factor(metadata$treatment)

# Set random seed for reproducibility
set.seed(123)

# Perform UMAP dimensionality reduction
umap_result <- umap(gene_expression_data, n_neighbors = 5, min_dist = 0.5)

# Extract UMAP coordinates and combine with metadata
umap_df <- data.frame(
  X1 = umap_result$layout[, 1],  # UMAP component 1
  X2 = umap_result$layout[, 2],  # UMAP component 2
  metadata
)

# Plot using ggplot2
ggplot(umap_df, aes(x = X1, y = X2, color = treatment)) +
  geom_point(size = 3) +
  labs(title = "GSE152620_UMAP",
       x = "UMAP Component 1", y = "UMAP Component 2") +
  theme_minimal()

################################################## PCA ############################3
# 
# library(WGCNA)
# library(DESeq2)
# library(GEOquery)
# library(tidyverse)
# #library(CorLevelPlot)
# library(gridExtra)
# 
# data <- gene_expression_data
# #data1 <- read.csv("Count_Data_final.csv", header = TRUE, row.names = 1)
# data <- as.data.frame(data)
# #metadata <- read.csv("Meta_Data_final.csv", header = TRUE, row.names = 1) 
# metadata <- metadata
# # Quality Control: detect outlier genes
# gsg <- goodSamplesGenes(data)
# summary(gsg)
# 
# # Remove outlier genes
# data <- data[gsg$goodGenes == TRUE,]

# Detect outlier samples using hierarchical clustering
htree <- hclust(dist(gene_expression_data), method = "average")
######
#######
# Plot hierarchical clustering with colors based on metadata groups
plot(htree, labels = metadata$Sample.Name, main = "Hierarchical Clustering Dendrogram", col=metadata$treatment)
rect.hclust(htree, k = length(unique(metadata$group)), border = 2:5)

