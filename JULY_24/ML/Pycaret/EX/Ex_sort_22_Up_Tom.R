
setwd("E:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL/June_Run/GSE233233_QZ_006")
setwd("E:/DWCT/M_Analysis_Ajit_Sir/JULY_24/ML/Pycaret")
library(readr)
data <- read_csv("Normalized_Salt_vs_control.csv")
View(data)

genes22 <- read.csv("22_UP_RNA.csv")
entrez <- read.csv("gene_entrez_UP.csv")
ID <- entrez[entrez$To %in% genes22$UP,]

write.csv(ID,file = "22_UP_Entrez.csv")


ex22_up <- data[data$...1 %in% ID$From,]


write.csv(ex22_up, file="QZ_006_EX_22_UP.csv")

###############################

data <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL/June_Run/GSE233233_Zhongza_9/Normalized_Salt_vs_control.csv")


ex22_up <- data[data$...1 %in% ID$From,]

write.csv(ex22_up, file="Zhongza_9_EX_22_UP.csv")

###############################################



data <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE148353/June_Run/June_Run_Variant/Normalized_Salt stress_vs_control_June.csv")

ex22_up <- data[data$...1 %in% ID$From,]

write.csv(ex22_up, file="Micro_TOM_EX_22_UP.csv")


###############################################


data <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE148353/June_Run/June_Run_Micro_TOM_EX/Normalized_Salt stress_vs_control_June.csv")

ex22_up <- data[data$...1 %in% ID$From,]

write.csv(ex22_up, file="Micro_TOM_EX_EX_22_UP.csv")

##################################



data <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE152620/June_Run/Normalized_Salt_vs_control.csv")


genes22 <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/JULY_24/ML/Pycaret/22_UP_RNA.csv")
entrez <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/JULY_24/ML/Pycaret/GSE152620_gene_name_29X.csv")

ID <- entrez[entrez$edit %in% genes22$UP,]

write.csv(ID,file = "22_UP_GSE152620_gene_name_29X.csv")



ex22_up <- data[data$...1 %in% ID$main,]

write.csv(ex22_up, file="Boludo_29_X_EX_22_UP.csv")



############################### Dot TO NON DOT ##########################








genes22 <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/JULY_24/ML/Pycaret/22_UP_RNA.csv")
entrez <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/JULY_24/ML/Pycaret/Dot_to_Non_Dot_ALL.csv")

ID <- entrez[entrez$Edit %in% genes22$UP,]

write.csv(ID,file = "22_UP_Dot_With_Non_Tom.csv")


######################################################3 At least three cultiver ##############

library(readr)
Three_Least <- read.csv("At_least_3_datasets_RNA_UP_data_frame.csv")

library(readxl)

gene_list1 <- read.csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL/June_Run/GSE233233_QZ_006/Converted_ID_UP_GSE233233_QZ_006.csv")
gene_list2 <- read.csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL/June_Run/GSE233233_Zhongza_9/Converted_ID_UP_GSE233233_Zhongza_9.csv")
gene_list3 <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE148353/June_Run/June_Run_Micro_TOM_EX/Converted_UP_Micro_TOM_EX.csv")
gene_list4 <- read_excel("E:/DWCT/M_Analysis_Ajit_Sir/GSE148353/June_Run/June_Run_Variant/Converted_UP_Micro_TOM.xlsx")
 

# Merge gene lists based on the column "To"
merged_gene_list <- merge(gene_list1, gene_list2, by = "To", all = FALSE)  # Inner join gene_list1 and gene_list2
merged_gene_list <- merge(merged_gene_list, gene_list3, by = "To", all = FALSE)  # Inner join with gene_list3
merged_gene_list <- merge(merged_gene_list, gene_list4, by = "To", all = FALSE)  # Inner join with gene_list4



# Read in the CSV files
gene_list1 <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL/June_Run/GSE233233_QZ_006/Converted_ID_UP_GSE233233_QZ_006.csv")[, "To", drop = FALSE]
gene_list2 <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE233233/DGE_FINAL/June_Run/GSE233233_Zhongza_9/Converted_ID_UP_GSE233233_Zhongza_9.csv")[, "To", drop = FALSE]
gene_list3 <- read_csv("E:/DWCT/M_Analysis_Ajit_Sir/GSE148353/June_Run/June_Run_Micro_TOM_EX/Converted_UP_Micro_TOM_EX.csv")[, "To", drop = FALSE]

# Read in the Excel file and select only the "To" column
gene_list4 <- read_excel("E:/DWCT/M_Analysis_Ajit_Sir/GSE148353/June_Run/June_Run_Variant/Converted_UP_Micro_TOM.xlsx")[, "To", drop = FALSE]

# List of data frames to merge
gene_lists <- list(gene_list1, gene_list2, gene_list3, gene_list4)

# Merge data frames based on the column "To"
merged_gene_list <- Reduce(function(x, y) merge(x, y, by = "To", all = FALSE), gene_lists)

# merged_gene_list now contains the merged data frame based on the "To" column

combined_gene_list <- rbind(gene_list1, gene_list2, gene_list3, gene_list4)


write.csv(combined_gene_list,file="PPI_At_Least_Three_UP.csv")

getwd()


list <- read_csv("PPI_At_Least_Three_UP.csv")
 

Three_Least <- read.csv("At_least_3_datasets_RNA_UP_data_frame.csv")


# Extract columns 2 to 17 of Three_Least where values in each column are present in list$V2
selected_genes <- lapply(Three_Least[, 2:17], function(col) col[col %in% list$Edit])

# Convert the selected genes into a data frame
selected_genes_df <- as.data.frame(selected_genes)


# Initialize a list to store the matched values for each column of Three_Least
matched_values <- list()

# Iterate over columns 2 to 17 of Three_Least
for (i in 2:17) {
  # Extract values from the first column of list that match with current column of Three_Least
  matched_values[[i]] <- list$Edit[list$Edit %in% Three_Least[[i]]]
}

# Convert matched_values into a data frame if needed
matched_values_df <- as.data.frame(matched_values)

# Print or use matched_values_df as needed



rr1 <- list[list$Edit %in% Three_Least$Zhongza_9..QZ_006..TOM_EX,]
rr1 <- rr1[-2]




