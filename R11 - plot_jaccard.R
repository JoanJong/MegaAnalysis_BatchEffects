
### 3/10/19
### Violin plot of Jaccard scores for all methods

#### Libraries and functions ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001")
source('all_functions.R')

#### Importing data ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D06 - bootstrap_jaccard")

### Import Jaccard scores for all methods
comb_jac_imp = readRDS('comb_jac_mat_Razbatch.rds')

### Import bootstrap matrix
boot_list = readRDS('boot_list_Razbatch.rds')

#### Violin plot of Jaccard scores for all methods ####

# Need to restructure the dataset first
# Right now, every column is a diff normalisation method
# But we want them to be in 1 column

# No need janitor if check.names = F in read.table!!

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D07 - plot_jaccard")

tidy_comb_jac = comb_jac_imp %>%
  as.data.frame() %>% 
  gather(key = method, value = jac_score)

tidy_comb_jac$method = as.factor(tidy_comb_jac$method)

n_methods = nlevels(tidy_comb_jac$method)
palette = colorRampPalette(c("darkgreen", "lightgreen"))

plot = ggplot(tidy_comb_jac, aes(x = reorder(method, -jac_score, FUN=median), y = jac_score))

plot + geom_violin(aes(fill = reorder(method, -jac_score, FUN=median)), size = 0.3) +
  #stat_summary(aes(x = method, y = jac_score), fun.data = 'mean_sdl', geom="errorbar", size = 0.7, color = 'red') + 
  stat_summary(fun.y="median", geom="point", size=1, color="black") +
  theme(legend.position="none", text = element_text(size = 10), 
        axis.text.x = element_text(angle = 90, hjust = 1, colour = 'black')) + 
  labs(x = 'Method', y = 'Jaccard Score') +
  scale_fill_manual(values = palette(n_methods))

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D07 - plot_jaccard")
ggsave('plot_jaccard_violin_Razbatch.jpeg')
