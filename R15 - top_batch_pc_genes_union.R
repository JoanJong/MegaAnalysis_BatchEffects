
#### 30/7/20
#### Get the union of all top 5% genes
#### Find the frequencies of each gene appearing
#### Select the genes with the highest freq

#### Dates modified: 
#### 7/8/20 - converted table to dataframe, got gene freq in desc order
#### 13/8/20 - added pair names which gene occurs for each gene
#### 20/1/21 - using combated raz megadata

## All functions and libraries
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001")
source('all_functions.R')

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D09 - batch_effect_pairwise_genes")
pair_gene_mat = readRDS('pair_gene_mat_Razbatch.rds')

# pair_gene_mat is alr top 5% genes

#### --- Get union of top batch PC genes --- ####

#### UPDATED AS OF 13/8/20

## Combine all top genes from the 12 pairings into 1 vector
## For each gene, specify the pairing it came from
## Output: 12636 x 2, cols are 'gene' and 'pair_name'

gene_union_df = pair_gene_mat %>% 
  as.vector() %>% 
  as.data.frame() %>% 
  mutate(pair_name = rep(colnames(pair_gene_mat), each = 1053)) %>% 
  set_colnames(c('gene', 'pair_name'))

## Want to find freq of genes occuring, while keeping the pairing they came from

gene_union_unique = unique(as.character(gene_union_df$gene)) # 6766 genes

gene_freq_pairing_df = cbind(gene_union_unique, as.data.frame(matrix(data = 0L, 
                                                                     nrow = length(gene_union_unique), 
                                                                     ncol = ncol(pair_gene_mat))),
                             stringsAsFactors = F)
colnames(gene_freq_pairing_df) = c('gene', colnames(pair_gene_mat))
rownames(gene_freq_pairing_df) = gene_union_unique
  
for (gene in gene_union_unique)
{
  ## Subset the gene union set for a certain gene
  gene_subset = gene_union_df[gene_union_df$gene == gene, ]
  
  ## Get all the pairings for a certain gene
  pairing_subset = as.character(gene_subset$pair_name)
  
  for (pair_name in pairing_subset)
  {
    gene_freq_pairing_df[gene, pair_name] = 1
  }
  
}

gene_freq_pairing_df %<>%
  add_column(rowSums(.[, -1]), .after = 'gene') %>% 
  set_colnames(c('gene', 'freq', colnames(pair_gene_mat)))

#### --- Export df --- ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D11 - top_batch_pc_genes_union")
saveRDS(gene_freq_pairing_df, 'gene_freq_pairing_df_Razbatch.rds')

#### --------------------------------------------------------------------- ####

## Find union among all gene sets for every pairwise dataset
## split -- split pair_gene_mat into list elements, where each element is a column form pair_gene_mat
## reduce -- combine list elements into single value by iteratively applying union_all function
## union_all -- union with no removal of duplicates
gene_union = Reduce(union_all, split(pair_gene_mat, col(pair_gene_mat)))


gene_union_freqdesc_df = gene_union %>% 
  table() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  set_colnames(c('gene', 'freq')) %>% 
  set_rownames(.$gene)
  
## Export gene frequencies
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D11 - top_batch_pc_genes_union")
saveRDS(gene_union_freqdesc_df, 'gene_union_freqdesc_df.rds')

## Combine gene frequencies and num of probes per gene
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D11 - top_batch_pc_genes_union")
gene_union_freqdesc_df = readRDS('gene_union_freqdesc_df.rds')

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D10 - probes_per_gene")
gene_num_of_probes_df = readRDS('gene_num_of_probes_df.rds')

gene_num_of_probes_df_subset = gene_num_of_probes_df %>% 
  filter(gene %in% rownames(gene_union_freqdesc_df)) 

gene_union_freq_probes_df = merge.data.frame(gene_union_freqdesc_df, 
                                  gene_num_of_probes_df_subset, 
                                  sort = F)

# colnames(gene_union_freq_probes_df) = c('gene',
#                                         'freq',
#                                         'raz_GSE40645',
#                                         'hubal_GSE83352',
#                                         'mercken_GSE87105',
#                                         'gonzalez_GSE98613')

## Export combined df
setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D11 - top_batch_pc_genes_union")
saveRDS(gene_union_freq_probes_df, 'gene_union_freq_probes_df.rds')



