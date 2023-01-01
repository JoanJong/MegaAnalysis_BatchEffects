## 20/6/19 AUTOMATED VERSION
## This R script compresses dataframes by collapsing rows with identical gene symbols with the mean expression value.
## It also does row-means imputation.

## ------- 1) IMPORT TXT FILES CONTAINING DATAFRAMES ------- ##

library(tools)
# Working directory should be one folder above /mapped_csv
setwd("C:/Users/User/OneDrive - Nanyang Technological University/Y2 S2/BS9001 Research Experience/R/win-library/aging_datasets")

# All files imported have col 1 as duplicating gene symbols

filepaths = list.files('./nn_mapped_txt', pattern = '*.txt', full.names = TRUE) # Get file path of dataset
filenames = list.files('./nn_mapped_txt', pattern = '*.txt', full.names = FALSE) # Get file name of dataset
filenames = lapply(filenames, file_path_sans_ext) # Just the filename without extension
all_df = lapply(filepaths, read.table, header=T, sep='\t') # Read all mapped dataframes
names(all_df) = filenames

print(any(is.na(all_df[[1]])))
print(any(is.na(all_df[[2]])))
print(any(is.na(all_df[[3]])))

counter = 1

## ------- 2) ROW MEAN IMPUTATION ------- ##

# Row means imputation
# Get index (row, col) of NA values
# For each row, find mean
# Insert mean into missing value

row_mean_imput = function(df) {
  rows_na = which(is.na(df), arr.ind=TRUE) # Find out which rows have NAs
  all_rows_na = unique(rows_na[,1]) # Get unique indices for rows with NAs
  rows_w_na = df[all_rows_na,] # Get rows with NAs
  
  row_means = rowMeans(rows_w_na[2:ncol(rows_w_na)], na.rm = TRUE) # Get row means
  
  for(i in 1:nrow(rows_w_na)){
    rows_w_na[i,is.na(rows_w_na[i,])] = row_means[i] # Substitute NAs with row mean
  }
  
  df[all_rows_na,] = rows_w_na  
  return(df)
}

## ------- 3) COMPRESSION ------- ##

for (df in all_df) {
  print(names(all_df[counter]))
  
  df = row_mean_imput(df)
  
  print(any(is.na(df)))
  View(df)
  
  mat_list = list()
  
  gene_symb_all = as.vector(df$Gene)  
  gene_symb = unique(gene_symb_all)
  
  for (gene in gene_symb) {
    df_subset = subset(df, Gene == gene, select = -c(1)) # Subset out only rows with corresponding gene symbol
    mat_subset = as.matrix(df_subset)                    # Convert dataframe to matrix
    mat_mean = apply(mat_subset, 2, mean)                # Find mean for each column
    mat_list[[gene]] = mat_mean                          # Add vector of means to mat_list
  }
    
  new_mat2 = do.call("rbind", mat_list)                  # Each gene symbol maps to 1 row of expression data
  
  file_name = paste(filenames[counter], '_cps.txt', sep = '')
  write.table(new_mat2, file=paste('./nn_compressed_df/', file_name, sep = ''), sep="\t", row.names=T, col.names=T, quote=F)
  
  print(any(is.na(new_mat2)))
  
  print(counter)
  counter = counter + 1
  }
