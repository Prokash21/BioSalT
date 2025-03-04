#
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
setwd("D:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL")

#load the count data

count_data <- read.csv("Count_Data_final.csv", header=TRUE,row.names = 1)
count_data<- count_data[,1:18]
colnames(count_data)
head(count_data)

#load the sample info
#sample_info1 <- read.csv("meta_keratinocyte.csv")

sample_info <- read.csv("Meta_Data_final.csv", header = TRUE,row.names = 1)
sample_info <- sample_info[1:18, ]

#colData1 <- read.csv("meta_keratinocyte.csv", header = T, sep = '\t', 
#                     stringsAsFactors = TRUE)
#colData <- read.csv("meta_keratinocyte.csv", header = TRUE,row.names = 1)
colnames(sample_info)
head(sample_info)
dim(sample_info)
#set factor levels
sample_info$agent <- factor(sample_info$agent)
sample_info$Time_point <- factor(sample_info$Time_point)
sample_info$tissue <- factor(sample_info$tissue)

# Convert non-integer values to integers in count data
count_data <- round(count_data)
head(count_data)

# Create a new count data object
new_count_data <- as.matrix(count_data)
head(new_count_data)


unique(sample_info$agent)
unique(sample_info$Time_point)
unique(sample_info$tissue)


dim(count_data)
dim(new_count_data)
dim(sample_info)

# Generate the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = new_count_data, colData = sample_info, design = ~ Time_point)

# Perform DESeq2 analysis
dds <- DESeq(dds)
head(dds)

#set the factor level
#dds$Treatment <- factor(dds$agent, levels = c ("control","Salt stress")) 


#filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]
dds

#set the referene for the treatment factor
dds$Time_point<- relevel(dds$Time_point , ref = "0h")
dds$Time_point

#perform the statistical tests to identify differentialy expressed genes
dds <- DESeq(dds)
head(dds)

#save the normalized counts
normalize_counts <- counts(dds,normalized=TRUE)
head(normalize_counts)
dim(normalize_counts)
write.csv(normalize_counts,"Salt stress_vs_0h.csv")



#Identify available coefficient names
coeff_names <- resultsNames(dds)

#Print the coefficient names
print(coeff_names)

#[[1] "Intercept"            "Time_point_12h_vs_0h" "Time_point_24h_vs_0h"
#[4] "Time_point_48h_vs_0h" "Time_point_6h_vs_0h"  "Time_point_96h_vs_0h"
######################################## "Time_point_12h_vs_0h"  ################################################

resLFC <- lfcShrink(dds, coef ="Time_point_12h_vs_0h", type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)

XX<-read.csv("Gene ID.csv")
top <- XX$Gene.ID
top <- as.character(top)
Y<-resLFC <- resLFC[top, ]

write.csv(Y, file = "LFC_12h.csv")


######################################## "Time_point_24h_vs_0h" ################################################
resLFC <- lfcShrink(dds, coef ="Time_point_24h_vs_0h", type = "apeglm")
resLFC<- as.data.frame(resLFC)
#change resLFC to a dataframe
resLFC1 <- as.data.frame(resLFC)

XX<-read.csv("Gene ID.csv")
top <- XX$Gene.ID
print(top)
top <- as.character(top)
Y<- resLFC1[top, ]

write.csv(Y, file = "LFC_24h.csv")


######################################## "Time_point_48h_vs_0h" ################################################
resLFC <- lfcShrink(dds, coef ="Time_point_48h_vs_0h", type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)

XX<-read.csv("Gene ID.csv")
top <- XX$Gene.ID
top <- as.character(top)
Y<-resLFC <- resLFC[top, ]

write.csv(Y, file = "LFC_48h.csv")
######################################## "Time_point_6h_vs_0h" ################################################
resLFC <- lfcShrink(dds, coef ="Time_point_6h_vs_0h", type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)

XX<-read.csv("Gene ID.csv")
top <- XX$Gene.ID
top <- as.character(top)
Y<-resLFC <- resLFC[top, ]

write.csv(Y, file = "LFC_6h.csv")
######################################## "Time_point_96h_vs_0h" ################################################

resLFC <- lfcShrink(dds, coef ="Time_point_96h_vs_0h", type = "apeglm")

#change resLFC to a dataframe
resLFC <- as.data.frame(resLFC)

XX<-read.csv("Gene ID.csv")
top <- XX$Gene.ID
top <- as.character(top)
Y <- resLFC[top, ]

write.csv(Y, file = "LFC_96h.csv")
















