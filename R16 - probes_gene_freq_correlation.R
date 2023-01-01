
#### 13/8/20
#### Find the correlation between number of probes per gene with freq of gene appearing in gene union set.
#### Does the number of probes a gene have affect its contribution to batch effect?

#### --- 1) Import the df with probe num and gene freq --- ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D11 - top_batch_pc_genes_union")

gene_probe_freq_imp = readRDS('gene_union_freq_probes_df.rds')

#### --- 2) Find the average probe num for each gene across the 4 datasets --- ####

## Indep variable is number of probes
## Dep variable is gene freq in gene union set
## x-axis (number of probes) sorted low to high

gene_probe_freq_mean = gene_probe_freq_imp %>% 
  mutate(mean_probe = rowMeans(.[, 3:6])) %>% 
  arrange(mean_probe) %>% 
  mutate_at(vars(freq), as.factor)

levels(gene_probe_freq_mean$freq)

#### --- 3) Plot gene freq against probe num --- ####

ggplot(gene_probe_freq_mean, aes(x = freq, y = mean_probe)) +
  geom_boxplot() +
  ## Set the number of tick marks
  scale_y_continuous(breaks = seq(0, (max(gene_probe_freq_mean$mean_probe) + 2), 2)) +
  ## Customise axis labels
  labs(x = 'Freq of Gene in Union Set', y = 'Mean Number of Probes per Gene')

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D13 - probes_gene_freq_corelation")

ggsave2('plot_probenum_agst_genefreq.jpeg')

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D13 - probes_gene_freq_corelation")

saveRDS(gene_probe_freq_mean, 'gene_probe_freq_mean.rds')
