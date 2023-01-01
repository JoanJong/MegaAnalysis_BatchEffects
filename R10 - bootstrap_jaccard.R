
### 25/10/19
# • Import normalized datasets *
# • Do bootstrapping to get samples -- 1000 rows of young vs old samples *
# • T-test between genes (rows) in young and old samples for each bootstrap *
# • Binarize p-values *
# • Export binarized p-values from T-test *
# • Pairwise comparison of binarized p-values between bootstraps (columns) and get Jaccard scores for each comparison (SX) *
# • Export Jaccard scores
# • Plot combined violin plot of Jaccard scores for all datasets and export
# • Sum up binary values for each gene (row) to get rowsums
# • Density line plot of rowsums for all datasets

#### Libraries and functions ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001")
source('all_functions.R')

#### 1) Data import and initialisation ####

#### 1/1/21 - Using exprs data with Raz already internally ComBated

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D08 - making_eset")
ori_eset_Razbatch = readRDS('ori_eset_Razbatch.rds')

class_fac = readRDS('class_factor.rds')
batch_fac = readRDS('batch_factor.rds')
gender_fac = pData(ori_eset)$sex

neat_dataset_names = readRDS('neat_BERM_names.rds')

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D04 - norm_combat_data")
normed_data = readRDS('normed_data_Razbatch.rds')

# Set number of bootstrap samples
num_of_boot = 1000

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001/different_age_thresholds/50")

df_import = read.table('dataset_5050.txt', sep='\t', header=T, stringsAsFactors = FALSE)

class_fac = as.factor(df_import[5, ]) # 66 old, 44 young
batch_fac = as.factor(df_import[2, ]) # 4 batches
gender_fac = as.factor(df_import[4, ]) # male, female

#### 2) Sampling columns ####

# Sample the columns from df_import

# Get GSM names of samples from young
df_imp_y_samples = df_import[1, class_fac == 'young']
# Get GSM names of samples from old
df_imp_o_samples = df_import[1, class_fac == 'old']

# Remove colnames
colnames(df_imp_y_samples) = NULL; colnames(df_imp_o_samples) = NULL
rownames(df_imp_y_samples) = NULL; rownames(df_imp_o_samples) = NULL

# Number of samples in each class
num_y_samples = length(df_imp_y_samples)
num_o_samples = length(df_imp_o_samples)

# Initialise matrix
sample_col_mat = matrix(NA, nrow = num_of_boot, ncol = ncol(df_import))

set.seed(123)

for (resample_idx in 1:num_of_boot)
{
  # Resample from young samples
  sample_col_mat[resample_idx, 1:num_y_samples] = unlist(sample(x = df_imp_y_samples, size = num_y_samples, replace = T))
  # Resample from old samples
  sample_col_mat[resample_idx, (num_y_samples+1):ncol(df_import)] = unlist(sample(x = df_imp_o_samples, size = num_o_samples, replace = T))
  
}

# Check if all samples are resampled at least once
setequal(as.character(df_import[1, ]), as.data.frame(table(sample_col_mat), stringsAsFactors = F)[, 1])

# Get gene names - 21053 genes
gene_names = rownames(df_import)[7:nrow(df_import)]

#### 3) Bootstrapping + T-test ####

# setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001/bootstrap_matrices")

boot_list = vector('list', length(dataset_names))
names(boot_list) = dataset_names

##### LONG TIME!!! #####
for (data_name in dataset_names)
{
  print(data_name)
  
  boot_list[[data_name]] = boot_t_function_v2(normed_data[[data_name]], sample_col_mat, num_of_boot, class_fac,
                                              num_y_samples, num_o_samples)
  
  # write.table(boot_list[[data_name]], file = paste0(data_name, '_tmat.txt'), sep="\t", col.names=T, row.names=T)
  
}

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D06 - bootstrap_jaccard")
saveRDS(boot_list, 'boot_list_Razbatch.rds')

#### Calculate Jaccard score ####

jac_list = vector('list', length(dataset_names))
names(jac_list) = dataset_names

# setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001/jaccard_scores")

for (data_name in dataset_names)
{
  jac_list[[data_name]] = jaccard_compare_func_v2(boot_list[[data_name]], num_of_boot = num_of_boot)
  
  # write.table(jac_list[[data_name]], file = paste0(data_name, '_jac.txt'), sep="\t", col.names=T, row.names=T)
  
}

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D06 - bootstrap_jaccard")
#saveRDS(jac_list, 'jac_list_Razbatch.rds')

#### Combine Jaccard scores into a dataframe ####

# Consolidate all jaccard scores for all selected methods

for(data_name in dataset_names)
{
  jac_vec = na.omit(as.vector(jac_list[[data_name]]))
  
  if (data_name == 'ori')
  {
    comb_jac_mat = jac_vec
  }
  
  else
  {
    comb_jac_mat = cbind(comb_jac_mat, jac_vec)
  }
  
}

colnames(comb_jac_mat) = neat_dataset_names

# Export
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D06 - bootstrap_jaccard")

saveRDS(comb_jac_mat, 'comb_jac_mat_Razbatch.rds')
