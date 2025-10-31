# Install and load the necessary packages
#install.packages("ggplot2")
#install.packages("ggVennDiagram")
library(ggplot2)
library(ggVennDiagram)
setwd("E:/DWCT/M_Analysis_Ajit_Sir/JULY_24/ggVenn")

data <-read.csv("UP_RNA_ALL.csv")

# Example gene lists
gene_list1 <- data$GSE233233_Zhongza_9
gene_list2 <- data$GSE233233_QZ_006
gene_list3 <- data$Micro_TOM_EX
gene_list4 <- data$Micro_TOM
gene_list5 <- data$X29_X_Boludo

# Create a list of gene lists
gene_lists <- list(
  Zhongza_9 = gene_list1,
  QZ_006 = gene_list2,
  TOM_EX = gene_list3,
  Micro_TOM = gene_list4,
  X29_X_Boludo = gene_list5
)

# Generate the Venn diagram
ggVennDiagram(gene_lists, label = "count")



########################################################################################
# Load the necessary packages
#install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# Example gene lists from your data
gene_list1 <- data$GSE233233_Zhongza_9
gene_list2 <- data$GSE233233_QZ_006
gene_list3 <- data$Micro_TOM_EX
gene_list4 <- data$Micro_TOM
gene_list5 <- data$X29_X_Boludo

# Create a list of gene lists
gene_lists <- list(
  Zhongza_9 = gene_list1,
  QZ_006 = gene_list2,
  TOM_EX = gene_list3,
  Micro_TOM = gene_list4,
  X29_X_Boludo = gene_list5
)

# Generate the Venn diagram
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("Zhongza_9", "QZ_006", "TOM_EX", "Micro_TOM", "X29_X_Boludo"),
  filename = NULL,  # to plot in RStudio
  output = TRUE,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = "black",
  fill = c("dodgerblue", "goldenrod1", "seagreen3", "darkorchid1", "red"),
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135, -135, 180),
  cat.dist = c(0.055, 0.055, 0.085, 0.085, 0.085),
  cat.fontfamily = "sans",
  rotation.degree = 1,
  euler.d = TRUE,
  scaled = TRUE
)

# Plot the Venn diagram
grid.draw(venn.plot)
X<- as.data.frame(venn.plot)



# Calculate the lengths of each gene list and intersections
lengths <- sapply(gene_lists, length)
combinations <- expand.grid(rep(list(c(FALSE, TRUE)), 5))

intersection_lengths <- apply(combinations, 1, function(row) {
  selected_lists <- gene_lists[row]
  if (all(!row)) {
    return(0)
  }
  Reduce(intersect, selected_lists)
})

# Create a data frame of lengths and intersection lengths
intersection_lengths_df <- data.frame(
  Combination = apply(combinations, 1, function(x) paste(names(gene_lists)[x], collapse = ", ")),
  Length = sapply(intersection_lengths, length)
)

# Display the data frame
print(intersection_lengths_df)


########################################

# Calculate the intersections
combinations <- expand.grid(rep(list(c(FALSE, TRUE)), length(gene_lists)))

# Function to find the intersection of selected gene lists
find_intersections <- function(row) {
  selected_lists <- gene_lists[row]
  if (sum(row) < 3) {  # Filter for intersections of at least 3 gene lists
    return(character(0))
  }
  Reduce(intersect, selected_lists)
}

# Apply the function to each combination
intersection_results <- apply(combinations, 1, find_intersections)

# Create a data frame for combinations and intersections
intersection_df <- data.frame(
  Combination = apply(combinations, 1, function(x) paste(names(gene_lists)[x], collapse = ", ")),
  Genes = sapply(intersection_results, function(x) paste(x, collapse = ", ")),
  stringsAsFactors = FALSE
)

# Filter out empty intersections and those with less than 3 sets
intersection_df <- intersection_df[intersection_df$Genes != "" & rowSums(combinations) >= 3, ]

# Display the data frame with intersections of at least 3 gene sets
print(intersection_df)

write.csv(intersection_df,file="At_least_3_datasets_RNA_UP.csv")

# Create a data frame for combinations and intersections
intersection_df <- data.frame(
  Combination = apply(combinations, 1, function(x) paste(names(gene_lists)[x], collapse = ", ")),
  Genes = sapply(intersection_results, function(x) paste(x, collapse = ", ")),
  stringsAsFactors = FALSE
)

# Filter out empty intersections and those with less than 3 sets
intersection_df <- intersection_df[intersection_df$Genes != "" & rowSums(combinations) >= 3, ]

# Split genes into a list
intersection_list <- strsplit(intersection_df$Genes, ", ")

# Find the maximum length of gene lists
max_length <- max(sapply(intersection_list, length))

# Create a matrix to hold the genes with NA padding
gene_matrix <- do.call(cbind, lapply(intersection_list, function(x) {
  c(x, rep(NA, max_length - length(x)))
}))

# Convert the matrix to a data frame and set column names
gene_df <- as.data.frame(gene_matrix, stringsAsFactors = FALSE)
colnames(gene_df) <- intersection_df$Combination
#colnames(dim(gene_df)) <- intersection_df$Combination

# Display the data frame
print(gene_df)

write.csv(gene_df,file="At_least_3_datasets_RNA_UP_data_frame.csv")



