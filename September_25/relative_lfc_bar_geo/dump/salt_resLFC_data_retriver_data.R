



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
gene_list <- c("Solyc02g084850.3.1", "Solyc06g067980.3.1",
               "Solyc06g067980.3.1"
               # , "URS00005CF03F_Y3",
               # "URS000032B6B6_U1", "URS00006CB9B5_U3",
               # "URS0000635FD4_U5", "URS00006F89B3_U8",
               # "URS000057BD5C_U35"
               )



###################

# take all annotated_resLFC files
csv_files <- c(
  "Buludo_29_X_resLFC_p_cut.csv"
  # "resLFC_Treatment_eHeLa_vs_eADMSC.csv",
  # "resLFC_Treatment_eMDAMB231_vs_eADMSC.csv",
  # "resLFC_Treatment_eA549_vs_eADMSC.csv",
  # "resLFC_Treatment_eH1975_vs_eADMSC.csv"
)

# # take all annotated_resLFC files
# csv_files <- c(
#   "resLFC_Treatment_eMCF7_vs_eBMMSC.csv",
#   "resLFC_Treatment_eHeLa_vs_eBMMSC.csv",
#   "resLFC_Treatment_eMDAMB231_vs_eBMMSC.csv",
#   "resLFC_Treatment_eA549_vs_eBMMSC.csv",
#   "resLFC_Treatment_eH1975_vs_eBMMSC.csv"
# )
###################

# condition names for ev 
column_names <- c("Buludo")
                  #, "eHeLa", "eMDAMB231", "eA549", "eH1975")

###################

file_name <- "Bulodo"

#control <- "eBMMSC"

# Call the function
result <- extract_all_metrics_multiple(gene_list, csv_files, column_names)

# Save the output
write.csv(result, file = paste0(file_name, "_resLFC_data_for_bar.csv"), row.names = FALSE)


