



######################

extract_all_metrics_multiple <- function(gene_list, csv_files, column_names) {
  # Initialize the data frame with gene names as the first column
  gene_frame <- data.frame(Gene = gene_list)
  
  # Loop over the provided CSV files and column names
  for (i in seq_along(csv_files)) {
    csv_file <- csv_files[i]
    column_prefix <- column_names[i]  # Prefix for condition-specific columns
    
    # Check if the file exists
    if (!file.exists(csv_file)) {
      cat("File", csv_file, "does not exist. Skipping...\n")
      next
    }
    
    # Read the CSV file
    data <- read.csv(csv_file)
    
    # Define the columns to extract
    selected_columns <- c("X", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
    
    # Filter the data for selected genes
    gene_data <- data[data$X %in% gene_list, selected_columns]
    
    # Match extracted values with gene list order
    matched_data <- t(sapply(gene_list, function(gene) {
      row <- gene_data[gene_data$X == gene, ]
      if (nrow(row) > 0) {
        as.numeric(row[1, -1])  # Exclude X column
      } else {
        rep(NA, length(selected_columns) - 1)  # Return NA if gene not found
      }
    }))
    
    # Convert to data frame
    matched_df <- as.data.frame(matched_data)
    colnames(matched_df) <- paste(column_prefix, c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"), sep = "_")
    
    # Add to main gene_frame
    gene_frame <- cbind(gene_frame, matched_df)
  }
  
  return(gene_frame)
}




################################################################################



# enter the gene list
gene_list <- c("URS0000007D24_Y4", "URS0000103047_Y1",
               "URS00004A2461_Y5", "URS00005CF03F_Y3",
               "URS000032B6B6_U1", "URS00006CB9B5_U3",
               "URS0000635FD4_U5", "URS00006F89B3_U8",
               "URS000057BD5C_U35")



###################

# take all annotated_resLFC files
csv_files <- c(
  "resLFC_Treatment_eMCF7_vs_eADMSC.csv",
  "resLFC_Treatment_eHeLa_vs_eADMSC.csv",
  "resLFC_Treatment_eMDAMB231_vs_eADMSC.csv",
  "resLFC_Treatment_eA549_vs_eADMSC.csv",
  "resLFC_Treatment_eH1975_vs_eADMSC.csv"
)

# take all annotated_resLFC files
csv_files <- c(
  "resLFC_Treatment_eMCF7_vs_eBMMSC.csv",
  "resLFC_Treatment_eHeLa_vs_eBMMSC.csv",
  "resLFC_Treatment_eMDAMB231_vs_eBMMSC.csv",
  "resLFC_Treatment_eA549_vs_eBMMSC.csv",
  "resLFC_Treatment_eH1975_vs_eBMMSC.csv"
)
###################

# condition names for ev 
column_names <- c("eMCF7", "eHeLa", "eMDAMB231", "eA549", "eH1975")

###################

control <- "eADMSC"
control <- "eBMMSC"

# Call the function
result <- extract_all_metrics_multiple(gene_list, csv_files, column_names)

# Save the output
write.csv(result, file = paste0(control, "_resLFC_data_for_heatmap_bar.csv"), row.names = FALSE)


################################################################################
# Heatmap generation
################################################################################

heatmap_data <- result[, grepl("log2FoldChange", colnames(result))]
row.names(heatmap_data) <- result$Gene
colnames(heatmap_data) <- column_names
heatmap_data[is.na(heatmap_data)] <- 0
heatmap_data[is.infinite(heatmap_data)] <- 0

# Install necessary packages if not already installed
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

# Load necessary libraries
library(pheatmap)
library(RColorBrewer)


# Function to create a heatmap
create_heatmap <- function(heatmap_data, save_path = "eBMMSC_heatmap_output") {
  # Convert to matrix
  heatmap_data <- as.matrix(heatmap_data)
  
  # Define color palette (from blue to red)
  color_palette <- colorRampPalette(c("#804a7e", "#f6f5fa", "#f29c5a"))(1000)
  
  # Define breaks for coloring
  max_abs_value <- max(abs(heatmap_data), na.rm = TRUE)
  breaks <- seq(-max_abs_value, max_abs_value, length.out = 1000)
  
  # Create heatmap
  heatmap_plot <- pheatmap(heatmap_data,
                           col = color_palette,
                           border_color = "#696880",
                           breaks = breaks,
                           clustering_distance_rows = "euclidean",
                           clustering_distance_cols = "euclidean",
                           clustering_method = "complete",
                           cellwidth = 24,
                           cellheight = 16,
                           fontsize_row = 8,
                           fontsize_col = 8,
                           angle_col = 45)
  
  # Modify column names to only show condition names
  colnames(heatmap_data) <- column_names
  
  # Save the heatmap in different formats
  pdf(file = paste0(save_path, ".pdf"), width = 10, height = 8)
  print(heatmap_plot)
  dev.off()
  
  png(file = paste0(save_path, ".png"), width = 800, height = 600)
  print(heatmap_plot)
  dev.off()
  
  tiff(file = paste0(save_path, ".tiff"), width = 800, height = 600, res = 300)
  print(heatmap_plot)
  dev.off()
}

# Generate heatmap
create_heatmap(heatmap_data)

