
### 15/10/19
### Outputs normalised datasets
### - Import dataset_5050.txt
### - Get batch and class factors from dataset_5050.txt
### - Conduct normalization + BECAs on dataset
### - Boxplots for datasets to see distribution

#### Libraries and functions ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001")
source('all_functions.R')

#### Data import and initialisation #####

#### Using exprs data with Raz already internally ComBated

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D08 - making_eset")
ori_eset_Razbatch = readRDS('ori_eset_Razbatch.rds')

class_fac = readRDS('class_factor.rds')
batch_fac = readRDS('batch_factor.rds')
gender_fac = pData(ori_eset)$sex

df_exp = exprs(ori_eset_Razbatch)

# Initialise list to contain all the datasets
dataset_names = c('ori', 'qall', 'qall_comall', 'qclass', 'qclass_comall',
                  'comall', 'znorm', 'zcom')
normed_data = vector('list', length(dataset_names))
names(normed_data) = dataset_names

#### Normalisation + BECA application #####
##### 1) - Original data w/o normalisation -- real_ori

normed_data$ori = df_exp

##### 2) - QN on all -- real_qall 

#real_log2_mat = as.matrix(log2(df_exp))
real_log2_mat = df_exp
normed_data$qall = real_log2_mat

# Do QN on whole dataset
normed_data$qall = normalize.quantiles(normed_data$qall)

rownames(normed_data[['qall']]) = rownames(df_exp)
colnames(normed_data[['qall']]) = colnames(df_exp)

##### 3) - QN on all, ComBat for all -- real_qall_comall 

normed_data$qall_comall = real_log2_mat

# Do QN on whole dataset
normed_data$qall_comall = normalize.quantiles(normed_data$qall_comall)

# Do ComBat on whole dataset
normed_data$qall_comall = ComBat(normed_data$qall_comall, batch_fac, par.prior = T, prior.plots = F)

rownames(normed_data[['qall_comall']]) = rownames(df_exp)
colnames(normed_data[['qall_comall']]) = colnames(df_exp)

##### 4) - Class-specific QN -- real_qclass 

normed_data$qclass = real_log2_mat

# Do QN by class
normed_data$qclass[, which(class_fac == 'old')] = normalize.quantiles(normed_data$qclass[, which(class_fac == 'old')])
normed_data$qclass[, which(class_fac == 'young')] = normalize.quantiles(normed_data$qclass[, which(class_fac == 'young')])

##### 5) - Class-specific QN, ComBat for all -- real_qclass_comall 

normed_data$qclass_comall = real_log2_mat

# Do QN by class
normed_data$qclass_comall[, which(class_fac == 'old')] = normalize.quantiles(normed_data$qclass_comall[, which(class_fac == 'old')])
normed_data$qclass_comall[, which(class_fac == 'young')] = normalize.quantiles(normed_data$qclass_comall[, which(class_fac == 'young')])

# Combat on all
normed_data$qclass_comall = ComBat(normed_data$qclass_comall, batch_fac, par.prior = T, prior.plots = F)

##### 6) - ComBat for all -- real_comall 

normed_data$comall = real_log2_mat

normed_data$comall = ComBat(normed_data$comall, batch_fac, par.prior = T, prior.plots = F)

##### 7) - Z-Normalization -- real_znorm 

normed_data$znorm = z_norm_mega(real_log2_mat)

##### 8) - Z-Normalisation, ComBat for all -- real_zcom

normed_data$zcom = real_log2_mat

normed_data$zcom = z_norm_mega(normed_data$zcom)

normed_data$zcom = ComBat(normed_data$zcom, batch_fac, par.prior = T, prior.plots = F)

#### Saving normalised datasets ####

#### 1/1/21
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D04 - norm_combat_data")

saveRDS(normed_data, 'normed_data_Razbatch.rds')