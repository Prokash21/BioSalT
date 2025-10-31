## ----pre,echo=FALSE,results='hide'--------------------------------------------
library(knitr)
opts_chunk$set(warning=FALSE,message=FALSE,cache=FALSE)

## ----loadLibrary--------------------------------------------------------------
library(GEOquery)


## -----------------------------------------------------------------------------
# If you have network access, the more typical way to do this
# would be to use this:
gds <- getGEO("GSE233233")

# Install GEOquery package if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

# Load the GEOquery package
library(GEOquery)

# Download and load the GEO dataset
gds <- getGEO("GSE217631", GSEMatrix = TRUE)


# Check the structure of the loaded data
str(gds)

# Extract the first dataset (if there are multiple)
gse148353 <- gds[[1]]
X <- as.data.frame(gse148353)

# View the metadata
metadata <- pData(gse148353)
head(metadata)


# Install GEOquery package if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

# Load GEOquery package
library(GEOquery)

# Download and load the GEO dataset
gds <- getGEO("GSE148353", GSEMatrix = TRUE)

# Extract the expression set
expression_set <- gds[[1]]

# Get the count matrix
count_matrix <- exprs(expression_set)

# Display the first few rows of the count matrix
head(count_matrix)

