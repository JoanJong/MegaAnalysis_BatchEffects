
#### 25/8/20
#### Determining which batch effect removal method gives the highest reproducibility in terms of
#### similarity of top significant genes.
#### 1) Pair datasets in different combinations
#### 2) For each dataset pair, apply the 7 batch effect removal methods
#### 3) Get top significant genes from each pair with each method
#### 4) For each method, compare top significant genes between dataset pairs -- get similarity

#### 22/1/21
#### Using ComBated-raz dataset

########## 1) Data import and initialisation ##########

## Import log2 non-normalised dataset
setwd("C:/Users/sorae/OneDrive - Nanyang Technological University/Y4 S1 (OFYP)/Codes/D08 - making_eset")# ori_eset = readRDS('ori_50_eset.rds') 
ori_eset = readRDS('ori_eset_Razbatch.rds') 

class_fac = as.factor(ori_eset@phenoData@data$age) # 66 old, 44 young
batch_fac = as.factor(ori_eset@phenoData@data$batch_id) # 4 batches
names(batch_fac) = rownames(ori_eset@phenoData@data)

setwd("C:/Users/sorae/OneDrive - Nanyang Technological University/Y4 S1 (OFYP)/Codes/D08 - making_eset")# ori_eset = readRDS('ori_50_eset.rds') 
unknown_factor = readRDS('unknown_fac_Razbatch.rds')

########## 2) Find pairings of dataset ##########

## 6 pairings of datasets
dataset_pairings = combinations(v = levels(batch_fac),
                                n = 4, r = 2, 
                                repeats.allowed = F) %>%
  apply(., 1, function(x) {paste(x[1], x[2], sep = '_')})

## Create list:
## Each list element is one dataset pair, e.g. gonzalez_hubal
## Within each list element, 7 different batch correction mtds used

paired_datasets = vector('list', length = length(dataset_pairings))
names(paired_datasets) = dataset_pairings

pairwiseNormBatchCorrectionMod2021 = function(input_eset, pair_name)
{
  
  ## Get exprs data from eset
  exprs_data = exprs(input_eset)
  
  pheno_data = pData(input_eset)
  
  ## Get individual batch names from the pairing name
  pair_name_split = strsplit(x = pair_name, split = '_')
  
  ## Subset the exprs data for the batches desired
  pair_exprs = exprs_data[, input_eset@phenoData$batch_id == pair_name_split[[1]][1] | input_eset@phenoData$batch_id == pair_name_split[[1]][2]]
  
  ## Batch and class factors
  batch_factor = as.factor(pheno_data$batch_id[pheno_data$geo_accession %in% colnames(pair_exprs)]) %>% 
    droplevels()
  class_factor = as.factor(pheno_data$age[pheno_data$geo_accession %in% colnames(pair_exprs)])
  
  print(unknown_factor)
  
  ## Names of the batch effect removal mtds
  mtd_names = c('ori', 'qall', 'qall_comall', 'qclass', 'qclass_comall',
                'comall', 'znorm', 'zcom')
  
  ## Initialise list which will contain exprs data treated w the 8 mtds
  normed_data_list = vector('list', length = length(mtd_names))
  names(normed_data_list) = mtd_names
  
  ########## 1) Log2 only - ori ##########
  
  normed_data_list[['ori']] = pair_exprs
  print('1) ori: ok')
  
  ########## 2) QN on all - qall ##########
  
  normed_data_list[['qall']] = normalize.quantiles(pair_exprs)
  
  ## QN function removes row and colnames
  rownames(normed_data_list[['qall']]) = rownames(pair_exprs)
  colnames(normed_data_list[['qall']]) = colnames(pair_exprs)
  
  print('2) qall: ok')
  
  ########## 3) QN on all, Combat on all - qall_comall ##########
  
  normed_data_list[['qall_comall']] = ComBat(normed_data_list[['qall']], 
                                             batch_factor, par.prior = T, prior.plots = F)

  rownames(normed_data_list[['qall_comall']]) = rownames(pair_exprs)
  colnames(normed_data_list[['qall_comall']]) = colnames(pair_exprs)
  
  
  print('3) qall_comall: ok')
  
  ########## 4) Class-specific QN - qclass ##########
  
  normed_data_list[['qclass']] = pair_exprs
  
  normed_data_list[['qclass']][, which(class_factor == 'old')] = normalize.quantiles(normed_data_list[['qclass']][, which(class_factor == 'old')])
  normed_data_list[['qclass']][, which(class_factor == 'young')] = normalize.quantiles(normed_data_list[['qclass']][, which(class_factor == 'young')])
  
  ## QN function removes row and colnames
  rownames(normed_data_list[['qclass']]) = rownames(pair_exprs)
  colnames(normed_data_list[['qclass']]) = colnames(pair_exprs)
  
  print('4) qclass: ok')
  
  ########## 5) Class-specific QN, Combat on all - qclass_comall ##########
  
  normed_data_list[['qclass_comall']] = ComBat(normed_data_list[['qclass']],
                                               batch_factor, par.prior = T, prior.plots = F)
  print('5) qclass_comall: ok')
  
  ########## 6) Combat on all - comall ##########
  
  normed_data_list[['comall']] = ComBat(pair_exprs, 
                                        batch_factor, par.prior = T, prior.plots = F)
  print('6) comall: ok')
  
  ########## 7) ZN on all - znorm ##########
  
  normed_data_list[['znorm']] = z_norm_mega(pair_exprs)
  
  print('7) znorm: ok')
  
  ########## 8) ZN on all, Combat in all - zcom ##########
  
  normed_data_list[['zcom']] = ComBat(normed_data_list[['znorm']], 
                                      batch_factor, par.prior = T, prior.plots = F)
  print('8) zcom: ok')
  
  return(normed_data_list)
}

# MOD
for (pair_name_loop in dataset_pairings)
{
  print(pair_name_loop)
  paired_datasets[[pair_name_loop]] = pairwiseNormBatchCorrectionMod2021(ori_eset, pair_name_loop)
  
}


setwd("C:/Users/sorae/OneDrive - Nanyang Technological University/Y4 S1 (OFYP)/Codes/D17 - batch_mtds_pair_datasets")# saveRDS(paired_datasets, 'paired_datasets.rds')

saveRDS(paired_datasets, 'paired_datasets_Razbatch.rds')

saveRDS(batch_fac, 'batch_factor_all.rds')
saveRDS(class_fac, 'class_factor_all.rds')
