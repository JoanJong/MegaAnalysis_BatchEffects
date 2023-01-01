
### 24/10/19
# • Import dataset_5050.txt *
# • Get batch, class and gender factors from dataset_5050.txt *
# • Conduct PCA on normalized data *
# • Pairwise T-test within each PC to determine significance of difference when separating samples by class, batch and gender *
# • Binarize p-values *
# • Export normalized datasets
# • Export real PCA data (pca_data$x)
# • Export binarized p-value table (binlit)
# • Export PC loadings (pca_data$rotation)

#### Libraries and functions ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001")
source('all_functions.R')

#### Data import and initialisation #####

#### 1/1/21 - Using exprs data with Raz already internally ComBated

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D08 - making_eset")
ori_eset_Razbatch = readRDS('ori_eset_Razbatch.rds')

class_fac = readRDS('class_factor.rds')
batch_fac = readRDS('batch_factor.rds')
gender_fac = pData(ori_eset)$sex

# Import all normalised datasets
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D04 - norm_combat_data")
normed_data = readRDS('normed_data_Razbatch.rds')

#### Conduct PCA on normalised datasets ####

pc_num = 110

pca_data_list = vector('list', length(dataset_names))
names(pca_data_list) = dataset_names

for (data_name in dataset_names)
{
  ## Cannot use pca_data_list[[data_names[idx]]], not sure why!
  
  pca_data_list[[data_name]] = binary_generate(as.data.frame(t(normed_data[[data_name]])), pc_num, 
                                                        batchfac = batch_fac, classfac = class_fac, genderfac = gender_fac)
}

#### Export data ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D05 - pca_pwt_loads")

saveRDS(pca_data_list, 'pca_data_list_Razbatch.rds')