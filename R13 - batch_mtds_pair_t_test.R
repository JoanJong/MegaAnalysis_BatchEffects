
#### 28/8/20
#### Doing t-test for each paired dataset to get significant genes
#### For 1 mtd, see if sig genes are similar across all possible paired datasets

#### 7/5/21
#### Use p-val adjusting to limit number of genes reported sig

########## 1) Data import and initialisation ##########

setwd("C:/Users/sorae/OneDrive - Nanyang Technological University/Y4 S1 (OFYP)/Codes/D17 - batch_mtds_pair_datasets")
paired_datasets = readRDS('paired_datasets_Razbatch.rds')

batch_fac_all = readRDS('batch_factor_all.rds')
class_fac_all = readRDS('class_factor_all.rds')

mtds = names(paired_datasets$gonzalez_hubal)
pair_names = names(paired_datasets)

p_threshold = 0.05

########## 2) T-test ##########

## Initialise list for containing p-values for genes after doing t-test
paired_data_t_test = vector('list', length = length(paired_datasets))
names(paired_data_t_test) = names(paired_datasets)

## Do T-test
for (pair_name_loop in pair_names)
{
  print(pair_name_loop)
  
  paired_dataset_loop = paired_datasets[[pair_name_loop]]
  
  ## Get class factor for this dataset pairing
  class_fac_paired = class_fac_all[names(class_fac_all) %in% colnames(paired_datasets[[pair_name_loop]]$ori)]
  
  ## Go through the 8 mtds for this dataset pairing
  for (mtd in mtds)
  {
    ## For each gene, do T-test on samples separated by class
    paired_data_t_test[[pair_name_loop]][[mtd]] = paired_dataset_loop[[mtd]] %>% 
      apply(1, function(x) {t.test(x[class_fac_paired == 'old'], #paired_datasets[['mercken_raz']][['qclass_comall']][, class_fac_paired == 'old']
                                   x[class_fac_paired == 'young'])$p.value}) %>% 
      
      ## Arrange genes according to p-value (ascending)
      sort()
  }
}

## Find number of young and old samples for each pairing of datasets
for (pair_name_loop in pair_names)
{
  print(pair_name_loop)
  
  paired_dataset_loop = paired_datasets[[pair_name_loop]]
  
  class_fac_paired = class_fac_all[names(class_fac_all) %in% colnames(paired_datasets[[pair_name_loop]]$ori)]
  
  print(table(class_fac_paired))
  
}

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D18 - batch_mtds_pair_t_test")
saveRDS(paired_data_t_test, 'paired_data_p_vals_Razbatch.rds')

########## 3) Get top sig genes using p-values ##########

paired_data_siggenes = vector('list', length = length(paired_data_t_test))
names(paired_data_siggenes) = names(paired_data_t_test)

## Extract genes with p <= 0.05
for (pair_name_loop in pair_names)
{
  print(pair_name_loop)
  
  paired_siggenes_loop = paired_data_t_test[[pair_name_loop]]
  
  for (mtd in mtds)
  {
    print(mtd)
    
    ### P-val adjust - BH, bonferroni, holm, fdr
    pvals = p.adjust(paired_data_siggenes[[pair_name_loop]][[mtd]], method = "bonferroni")
    
    print(pvals)
    
    #paired_data_siggenes[[pair_name_loop]][[mtd]] = paired_siggenes_loop[[mtd]][paired_siggenes_loop[[mtd]] <= p_threshold] 
  }

}

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D18 - batch_mtds_pair_t_test")
saveRDS(paired_data_siggenes, 'paired_data_siggenes_Razbatch.rds')

########## 4) Calculate Jaccard similarity between sig genes across all dataset pairs for each mtd ##########

jaccardGeneMembership = function(gene_set_1, gene_set_2)
{
  jac_score = length(intersect(gene_set_1, gene_set_2)) / length(union(gene_set_1, gene_set_2))
  
  return(jac_score)
}

## For each mtd, find Jaccard score across all paired datasets
pairing_combi_vec = combinations(n = length(pair_names), 
                                 r = 2,
                                 v = pair_names,
                                 repeats.allowed = F) %>% 
  apply(1, function(x) {paste(x[1], x[2], sep = '&')})

paired_data_jac = as.data.frame(matrix(data = NA,
                                       nrow = length(mtds),
                                       ncol = length(pairing_combi_vec)))
colnames(paired_data_jac) = pairing_combi_vec
rownames(paired_data_jac) = mtds

for (mtd in mtds)
{
  
  print(mtd)
  
  for (pair_combi_loop in pairing_combi_vec)
  {
    
    print(pair_combi_loop)
    
    pair_name_loop = unlist(strsplit(pair_combi_loop, split = '&'))
    
    gene_set_1 = names(paired_data_siggenes[[pair_name_loop[1]]][[mtd]])
    gene_set_2 = names(paired_data_siggenes[[pair_name_loop[2]]][[mtd]])
    
    paired_data_jac[mtd, pair_combi_loop] = jaccardGeneMembership(gene_set_1 = gene_set_1,
                                                                  gene_set_2 = gene_set_2)
    
  }
  
}

paired_data_jac = t(paired_data_jac) %>% 
  as.data.frame() %>% 
  cbind(rownames(.), .)

colnames(paired_data_jac) = c('pairs', colnames(paired_data_jac)[-1])

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D18 - batch_mtds_pair_t_test")
saveRDS(paired_data_jac, 'paired_data_jac_Razbatch.rds')

########## 5) Boxplot for Jaccard score against mtds ##########

paired_data_jac_long = gather(data = paired_data_jac, 
                              key = batch_mtd,
                              value = jac_score,
                              ori:zcom,
                              factor_key=TRUE)

levels(paired_data_jac_long$batch_mtd) = c('None',
                                           'QN',
                                           'QN + Combat',
                                           'CSQN',
                                           'CSQN + Combat',
                                           'Combat',
                                           'ZN',
                                           'ZN + Combat')

ggplot(data = paired_data_jac_long, aes(x = batch_mtd, y = jac_score)) +
  geom_boxplot() +
  labs(x = 'Batch Effect Removal Method', 
       y = 'Jaccard Score', 
       title = 'Boxplot of Jaccard Scores for Batch Effect Removal Methods') +
  theme(axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        plot.title = element_text(hjust = 0.5))

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D18 - batch_mtds_pair_t_test")
ggsave('jaccard_paired_siggenes_Razbatch.jpg')

paired_data_medians = paired_data_jac_long %>% 
  group_by(batch_mtd) %>% 
  summarise(Median = median(jac_score), Std=sd(jac_score))

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D18 - batch_mtds_pair_t_test")
saveRDS(paired_data_medians, 'paired_data_medians_Razbatch.rds')

#### 11/9/20 Want to split boxplot into the different pairings

paired_data_jac_long %<>%
  ## Label for each mtd. Shld only show label on the end of the lines
  ## At other data points, is NA label
  mutate(label = if_else(pairs == 'hubal_raz&mercken_raz', as.character(batch_mtd), NA_character_))

# d_ends = paired_data_jac_long %>% 
#   group_by(batch_mtd) %>% 
#   top_n(1, jac_score) %>% 
#   pull(jac_score)

plot_aesthetics = theme(axis.text.x = element_text(angle = 90, 
                                                   vjust = 0.5, 
                                                   hjust = 1, 
                                                   size = 15),
                        axis.line.x = element_line(color="black", size = 0),
                        axis.line.y = element_line(color="black", size = 0),
                        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Batch Effect Removal Method',
       y = 'Jaccard Score',
       title = 'Similarity in Significant Gene Sets between Dataset Pairings')

## X-axis is dataset pairs name
ggplot(data = paired_data_jac_long, aes(x = pairs, y = jac_score, group = batch_mtd, colour = batch_mtd)) + 
  geom_line(size = 1.2) + 
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1, 
                                   size = 15),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = 'Batch Effect Removal Method',
       y = 'Jaccard Score',
       title = 'Similarity in Significant Gene Sets between Dataset Pairings')

## X-axis is dataset pairs name (to modify)
ggplot(data = paired_data_jac_long, aes(x = pairs, y = jac_score, group = batch_mtd, colour = batch_mtd)) + 
  geom_line(size = 1.2) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.line.x = element_line(color="black", size = 0),
        axis.line.y = element_line(color="black", size = 0),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.ticks = element_line(size = 1),
        axis.ticks.length.x = unit(1, 'mm')) +
  labs(x = NULL,
       y = NULL,
       title = 'Similarity in Significant Gene Sets between Dataset Pairings')

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D18 - batch_mtds_pair_t_test")
ggsave('jaccard_paired_data_siggenes_bypairs_empty_Razbatch.jpeg')

## X-axis is batch effect removal method name
ggplot(data = paired_data_jac_long, aes(x = batch_mtd, y = jac_score, group = pairs, colour = pairs)) + 
  geom_line(size = 1.2) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 15))

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D18 - batch_mtds_pair_t_test")
ggsave('jaccard_paired_data_siggenes.jpeg')
