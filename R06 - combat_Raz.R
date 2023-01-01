
#### 1/1/21 ####
#### Doing ComBat within Raz dataset first before combining with other datasets
#### to form mega-dataset.

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D08 - making_eset")

ori_eset = readRDS('ori_50_eset.rds')

class_fac = readRDS('class_factor.rds')
batch_fac = readRDS('batch_factor.rds')
gender_fac = pData(ori_eset)$sex

ori_exprs = exprs(ori_eset)

sample_names = rownames(pData(ori_eset))
raz_sample_names = names(batch_fac[batch_fac == 'raz'])

#### 1) Create Raz batch factor ####

group_1 = paste0('GSM998',
                 c(797, 798, 805, 806, 808, 810, 811, 813, 816,
                   802, 804, 809, 815, 820, 821, 822))


# Append group 1 with rest of Raz samples
unknown_fac_names = c(group_1, raz_sample_names[!raz_sample_names %in% group_1])

# Assign as 1 or 2 in the Raz batch factor
unknown_fac = factor(c(rep(1, length(group_1)), rep(2, length(raz_sample_names) - length(group_1))))

names(unknown_fac) = unknown_fac_names

# Reorder the Raz batch factor such that the sample names match that in the whole batch factor
unknown_fac = unknown_fac[order(match(names(unknown_fac), raz_sample_names))]

all(names(unknown_fac) == raz_sample_names)

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D08 - making_eset")
saveRDS(unknown_fac, 'unknown_fac_Razbatch.rds')

#### 2) Internally ComBat Raz ####

ori_exprs_mod = ori_exprs

ori_exprs_mod[, colnames(ori_exprs_mod) %in% raz_sample_names] = ComBat(ori_exprs_mod[, colnames(ori_exprs_mod) %in% raz_sample_names],
                                                                        unknown_fac, par.prior = T, prior.plots = F)

#### 3) Export modified ori_exprs ####

ori_eset_Razbatch = ori_eset

exprs(ori_eset_Razbatch) = ori_exprs_mod

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/D08 - making_eset")

# saveRDS(ori_eset_Razbatch, 'ori_eset_Razbatch.rds')

saveRDS(ori_eset_Razbatch, 'ori_eset_Razbatch.rds')


#### PCA scatterplot of Raz to figure out unknown factor ####

class_batch_fac = paste(class_fac, batch_fac)

pca_razbatch = prcomp(t(ori_exprs_mod[, colnames(ori_exprs_mod) %in% raz_sample_names]),
                      scale = F,
                      center = F)

autoplot(pca_razbatch)

ggplot(data = pca_razbatch$x[, c(1, 2)], 
       aes(x = PC1, y = PC2, label = rownames(pca_razbatch$x))) +
  geom_point(aes(shape = class_fac[names(class_fac) %in% raz_sample_names], 
                 colour = class_batch_fac[names(class_fac) %in% raz_sample_names],
                 fill = 'black'), size = 2.5) +
  geom_point(shape = class_fac[names(class_fac) %in% raz_sample_names],
             size = 2.5, colour = "black") +
  geom_label_repel() +
  # xlab(autoplot(pca_data_list$ori[[2]])$labels$x) +
  # ylab(autoplot(pca_data_list$ori[[2]])$labels$y) +
  scale_shape_manual(values = c(16, 17)) +
  scale_colour_manual(values = c('firebrick', 'orange')) +
  theme(legend.position = 'none')

###

pca_razbatch = prcomp(t(ori_exprs[, colnames(ori_exprs) %in% raz_sample_names]),
                      scale = F,
                      center = F)

ggplot(data = pca_razbatch$x[, c(1, 2)], 
       aes(x = PC1, y = PC2, label = rownames(pca_razbatch$x))) +
  geom_point(aes(shape = class_fac[names(class_fac) %in% raz_sample_names], 
                 colour = class_batch_fac[names(class_fac) %in% raz_sample_names],
                 fill = 'black'), size = 2.5) +
  geom_point(shape = class_fac[names(class_fac) %in% raz_sample_names],
             size = 2.5, colour = "black") +
  geom_label_repel() +
  xlab(autoplot(pca_razbatch)$labels$x) +
  ylab(autoplot(pca_razbatch)$labels$y) +
  scale_shape_manual(values = c(16, 17)) +
  scale_colour_manual(values = c('firebrick', 'orange')) +
  theme(legend.position = 'none')


