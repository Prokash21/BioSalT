
library(DESeq2)
library(pheatmap)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
getwd()
#set the working directory
setwd("D:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL")

#load the count data

count_data <- read.csv("Count_Data_final.csv", header=TRUE,row.names = 1)

colnames(count_data)
head(count_data)

#load the sample info
sample_info1 <- read.csv("Meta_Data_final_agent.csv")
sample_info <- read.csv("Meta_Data_final_agent.csv", header = TRUE,row.names = 1)
colData1 <- read.csv("Meta_Data_final_agent.csv", header = T, sep = '\t', 
                     stringsAsFactors = TRUE)
colData <- read.csv("Meta_Data_final_agent.csv", header = TRUE,row.names = 1)
colnames(sample_info)
head(sample_info)

#set factor levels
sample_info$agent <- factor(sample_info$agent)
sample_info$Cell_type <- factor(sample_info$Cell_type)

# Convert non-integer values to integers in count data
count_data <- round(count_data)
head(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)
head(new_count_data)


unique(sample_info$agent)
unique(sample_info$Cell_type)

all(rownames(count_data) == rownames(new_count_data))
all(colnames(count_data) == rownames(sample_info))
new_count_data <- new_count_data[rownames(sample_info), ]
X <- colnames(count_data)
Y <- rownames(sample_info)

X <- colnames(count_data)
Y <- rownames(sample_info)

# Find the names that are in X but not in Y
names_in_X_not_in_Y <- setdiff(X, Y)

# Find the names that are in Y but not in X
names_in_Y_not_in_X <- setdiff(Y, X)

# Print the results
print(names_in_X_not_in_Y)
print(names_in_Y_not_in_X)


# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~ agent)

# Perform DESeq2 analysis
dds <- DESeq(dds)
head(dds)

#set the factor level
dds$Treatment <- factor(dds$Treatment, levels = c ("mock","MPXV clade I infected", "MPXV clade IIb infected", "MPXV clade IIa infected")) 


#filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds

#set the referene for the treatment factor
dds$Treatment <- relevel(dds$Treatment , ref = "mock")
dds$Treatment

#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)

#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"normalized_counts_keratinocyte_I_vs_mock.csv")
