#### 22 Oct 21
##### Plot PCA scatterplots for each BERM-treated dataset

setwd("C:/Users/sorae/OneDrive - Nanyang Technological University/Y4 S1 (OFYP)/Codes/D05 - pca_pwt_loads")
pca_data_list = readRDS('pca_data_list_Razbatch.rds') # Raz-combat

setwd("C:/Users/sorae/OneDrive - Nanyang Technological University/Y4 S1 (OFYP)/Codes/D08 - making_eset")
ori_eset_Razbatch = readRDS('ori_eset_Razbatch.rds')

class_fac = readRDS('class_factor.rds')
batch_fac = readRDS('batch_factor.rds')
gender_fac = pData(ori_eset_Razbatch)$sex

#### Plot symbol ####
## 1 is circle (old), 3 is crosshair (young)

class_batch_fac = paste(class_fac, batch_fac)

neat_sample_labels = paste(rep(c('Old', 'Young'), each = 4), rep(c('Gonzalez', 'Hubal', 'Mercken', 'Raz'), times = 2))

colour_names = c('firebrick', 'orange',
                 'forestgreen', 'dodgerblue4',
                 'firebrick1', 'khaki',
                 'chartreuse2', 'deepskyblue')

setwd("C:/Users/sorae/OneDrive - Nanyang Technological University/Y4 S1 (OFYP)/Codes/D05 - pca_pwt_loads")

## Plot PCA

## Change 'ori' to another BERM name to generate PCA scatterplot for that BERM
## Also change plot title and file name accordingly
pcadata = pca_data_list$ori

plot_pca = ggplot(data = as.data.frame(pcadata[[2]]$x[, 1:2]), aes(x = PC1, y = PC2, label = rownames(pcadata[[2]]$x))) +
  geom_point(aes(shape = class_fac, colour = class_batch_fac, fill = 'black'), size = 2.5) +
  geom_point(shape = class_fac, size = 2.5, colour = "black") +
  # geom_label_repel() +
  scale_shape_manual(name = 'Legend',
                     labels = neat_sample_labels,
                     values = c(16, 17)) +
  scale_colour_manual(name = 'Legend',
                      labels = neat_sample_labels,
                      values = colour_names) +
  xlab(autoplot(pcadata[[2]])$labels$x) +
  ylab(autoplot(pcadata[[2]])$labels$y) +
  scale_x_continuous(breaks = seq(round(min(pcadata[[2]]$x[, 1]), digits = -1) - 10,
                                  round(max(pcadata[[2]]$x[, 1]), digits = -1) + 10, 5)) + # set PC1 range and interval
  scale_y_continuous(breaks = seq(round(min(pcadata[[2]]$x[, 2]), digits = -1) - 10,
                                  round(max(pcadata[[2]]$x[, 2]), digits = -1) + 10, 5)) + # set PC2 range and interval
  theme(legend.position = 'right') +
  guides(fill = FALSE,
         shape = FALSE,
         colour = guide_legend(override.aes = list(shape = rep(c(16, 17), each = 4)))) +
  ggtitle('Scatterplot of PC1 and PC2 for None-treated Dataset')

plot_pca

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D19 - plot_pca_scatterplot/ComBat-Raz")
ggsave('plot_ori_Razbatch.jpg')

dev.off()