
## 8/7/20

## 1) Import non-normalised dataset
## 2) For each class, get different combinations of batches
##    ○ Classes: young, old
##    ○ Batches: gon, mer, raz, hub
## 3) For each pairing, 
##    ○ Use binary_generate() function to do PCA and get binarised p-values and loadings
##    ○ Use get_genes_pc_loads() to get gene names for each PC
##    ○ Select only batch PCs
##    ○ Check tailed-ness of loading distribution -- if single-tailed, take top 5%, if double-tailed, take 2.5% from both ends
## 4) Find the intersection of genes between all pairings of the batches

#### 20/1/21
#### Using combated Raz data

#### --- 1) Import non-normalised dataset --- ####

## All functions and libraries
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001")
source('all_functions.R')

## Import non-normalised dataset
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D08 - making_eset")

# Combated-Raz non-norm megadata
ori_eset = readRDS('ori_eset_Razbatch.rds') 

#### --- 2) Different pairings of same-class-differnt-batch samples in non-normalised dataset --- ####

## Get combinations of batch names

batch_names = levels(pData(ori_eset)$'batch_id') # Get the names of the four datasets

pair_combi = combinations(4, 2, batch_names, repeats.allowed = F)

pair_combi = rbind(cbind(rep('young', 6), pair_combi), cbind(rep('old', 6), pair_combi))

#### --- 3) Generate PCA data --- ####

# 66 old, 44 young

## List to store PCA data for each pairwise combo

pair_names = apply(pair_combi, 1, function(x) {paste(x[1], x[2], x[3], sep = '_')})

pair_pcabin_list = vector('list', length(pair_names))
names(pair_pcabin_list) = pair_names

pair_loads_list = vector('list', length(pair_names))
names(pair_loads_list) = pair_names

## PCA
# Do PCA on each pairing
# Need customised batchfac, etc for each pairing

exprs_eset = exprs(ori_eset)

for (pair_name in pair_names)
{

  pair_name_split = strsplit(x = pair_name, split = '_')
  
  pair_exprs = exprs_eset[, ori_eset@phenoData$age == pair_name_split[[1]][1] & (ori_eset@phenoData$batch_id == pair_name_split[[1]][2] | ori_eset@phenoData$batch_id == pair_name_split[[1]][3])]
  
  pc_num = ncol(pair_exprs) # No of samples for one pair
  
  print(pc_num)
  
  # Need to get batch, class and gender factors unique for each pair
  # Have to do list then unlist because using c() will convert factors into integers
  
  factors_pair = pData(ori_eset)[ori_eset@phenoData$age == pair_name_split[[1]][1] & (ori_eset@phenoData$batch_id == pair_name_split[[1]][2] | ori_eset@phenoData$batch_id == pair_name_split[[1]][3]), ]
  
  batch_fac_pair = factors_pair$batch_id
  gender_fac_pair = factors_pair$sex
  
  ## Doing PCA
  
  binary_tab = data.frame('batch_bin' = rep(0, pc_num), 
                          'gender_bin' = rep(0, pc_num),
                          row.names = paste('PC', 1:pc_num, sep = ''))
  literal_tab = data.frame('batch_lit' = rep(0, pc_num),
                           'gender_lit' = rep(0, pc_num),
                           row.names = paste('PC', 1:pc_num, sep = ''))
  
  pca_data = prcomp(t(pair_exprs), scale = T, center = T) 
  pca_df = data.frame(pca_data$x[, 1:pc_num])
  
  ## Doing t-test for every PC
  
  for (pc in 1:pc_num) # Runs through each PC
  {
    
    # Can't use levels() because returns all 4 dataset names
    batch_t = t.test(pca_df[, pc][which(batch_fac_pair == pair_name_split[[1]][2])], 
                     pca_df[, pc][which(batch_fac_pair == pair_name_split[[1]][3])])
    
    gender_t = t.test(pca_df[, pc][which(gender_fac_pair == 'male')], 
                      pca_df[, pc][which(gender_fac_pair == 'female')])
  
    # Binarising p-vals
    # If p-val <= 0.05, assign as 1 (significant). Else, assign as 0 (not significant)
    
    if (batch_t$p.value <= 0.05)
    {
      binary_tab[pc, 'batch_bin'] = 1
    }
    
    if (gender_t$p.value <= 0.05)
    {
      binary_tab[pc, 'gender_bin'] = 1
    }
    
    # Insert the non-binarised p-vals into literal_tab matrix
    
    literal_tab[pc, 'batch_lit'] = batch_t$p.value
    literal_tab[pc, 'gender_lit'] = gender_t$p.value
    
  }
  
  pair_pcabin_list[[pair_name]] = cbind(binary_tab, literal_tab)
  
  pair_loads_list[[pair_name]] = pca_data
  
}

## Export PCA data and PCA t-test results

# setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/batch_effect_pairwise_genes")
# saveRDS(file = 'pair_loads_list.rds', pair_loads_list)
# saveRDS(file = 'pair_pcabin_list.rds', pair_pcabin_list)

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D09 - batch_effect_pairwise_genes")
saveRDS(file = 'pair_loads_list_Razbatch.rds', pair_loads_list)
saveRDS(file = 'pair_pcabin_list_Razbatch.rds', pair_pcabin_list)

#### --- 4) Get gene set for batch-correlated PC  --- ####

## Select first batch-correlated PC for each pair

first_batch_pc = rep(0, length(pair_names))

# first_batch_pc = paste0('PC', c(2, 2, 2, 1, 2, 2, 1, 2, 2, 1, 1, 2))

names(first_batch_pc) = pair_names

all(names(first_batch_pc) == names(pair_pcabin_list)) # TRUE


for (pair_name in pair_names)
{
  first_batch_pc[pair_name] = rownames(pair_pcabin_list[[pair_name]])[pair_pcabin_list[[pair_name]]$batch_bin == 1][1]
}

## Subset out loadings for the selected batch-correlated PC and get gene list

pair_gene_list = vector('list', length(pair_names))
names(pair_gene_list) = pair_names

for (pair_name in pair_names)
{
  
  print(pair_name)
  
  if (is.na(first_batch_pc[pair_name]) == T)
  {
    pair_gene_list[[pair_name]] = NA
    next()
  }
  
  ## 1) Get proportional PC loadings for each PC
  
  abs_loads = abs(pair_loads_list[[pair_name]]$rotation[, first_batch_pc[pair_name]]) # Convert all loadings to absolute values
  pc_loads_prop = abs_loads / sum(abs_loads) # Divide loadings by colSums to get proportional loadings
  
  ## 2) Gene names in decreasing order of rotation for each PC
  
  pair_gene_list[[pair_name]] = names(pc_loads_prop[order(pc_loads_prop, decreasing = T)])
  
  #rownames(pair_gene_list[[pair_name]]) = NULL
  
  # ## 3) Actual PC loadings arranged according to gene order
  # 
  # arr_load = loads
  # 
  # for (index in 1:ncol(pc_genes))
  # {
  #   # Gene vector for rows to be ordered by
  #   target = pc_genes[, index]
  #   
  #   # Reaarange rows for each column
  #   arr_load[, index] = arr_load[match(target, rownames(arr_load)), index]
  # }
  
}

## Export whole gene lists for all pairs

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D09 - batch_effect_pairwise_genes")
saveRDS(pair_gene_list, 'pair_gene_list_Razbatch.rds')

## Pick out top genes for each pair

top_batch_genes = round(nrow(exprs(ori_eset)) * 0.05)

pair_gene_mat = matrix(data = numeric(), nrow = top_batch_genes, ncol = length(pair_names))
colnames(pair_gene_mat) = pair_names

for (pair_name in pair_names)
{
  pair_gene_mat[, pair_name] = pair_gene_list[[pair_name]][1:top_batch_genes]
}

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D09 - batch_effect_pairwise_genes")

saveRDS(pair_gene_mat, 'pair_gene_mat_Razbatch.rds')

## Finding intersection of genes for all 12 pairings

Reduce(intersect, split(pair_gene_mat, col(pair_gene_mat))) # 31 genes


#### --- 5) Jaccard scoring --- ####

jaccard_similar_genes_func = function(x, y)
{
  jaccard_index = length(intersect(x, y)) / (length(x)+ length(y) - length(intersect(x, y)))
  
  if (is.nan(jaccard_index)) # Sometimes all values = 0 and jaccard index will give NaN
  {
    jaccard_index = NA
  }
  return(jaccard_index)
  
}

jaccard_mat = matrix(data = NA, nrow = length(pair_names), ncol = length(pair_names))
colnames(jaccard_mat) = pair_names
rownames(jaccard_mat) = colnames(jaccard_mat)
cols = c(1:length(pair_names))

for (row in 1:length(pair_names)) # Go through every row
{
  print(row)
  to_remove = c(1:row)
  cols = cols[! cols %in% to_remove]
  
  jaccard_mat[row, row] = 1
  
  for (col in cols) # Go through every pair for one row
  {
    jaccard_mat[row, col] = jaccard_similar_genes_func(pair_gene_mat[, row], pair_gene_mat[, col])
  }
  
}

jac_complete = na.omit(c(jaccard_mat))
jac_complete = jac_complete[!jac_complete == 1]

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/batch_effect_pairwise_genes")

saveRDS(first_batch_pc, 'first_batch_pc_Razbatch.rds')

saveRDS(pair_gene_mat, 'pair_gene_mat_Razbatch.rds')
