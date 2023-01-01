## 5/9/19
## Trying out different age thresholds to see if the sig genes reported are the same

library(GEOquery)

## Import dataset with metadata

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y2 S2/BS9001 Research Experience/R/win-library/aging_datasets")

sorted_comb = read.table('sorted_comb_data.txt', stringsAsFactors = F, sep = '\t', header = T)

View(sorted_comb)

## Add ages into the dataset

colnames(sorted_comb)

#### GSE40645 - actual age given ####

gse40645 = getGEO('GSE40645')
View(gse40645$GSE40645_series_matrix.txt.gz@phenoData@data)

gse40645_age_char = as.character(gse40645$GSE40645_series_matrix.txt.gz@phenoData@data[, 1])
gse40645_age_num_1 = sapply(gse40645_age_char[1:9], function(x){as.numeric(substring(x, first = 3, last = 4))})
gse40645_age_num_2 = sapply(gse40645_age_char[10:29], function(x){as.numeric(substring(x, first = 4, last = 5))})
gse40645_age = c(gse40645_age_num_1, gse40645_age_num_2)
names(gse40645_age) = rownames(gse40645$GSE40645_series_matrix.txt.gz@phenoData@data)

#### GSE83352 - actual age given ####
# filter out irrelevant samples

gse83352 = getGEO('GSE83352')
View(gse83352$GSE83352_series_matrix.txt.gz@phenoData@data)

# Get names of relevant samples
gse83352_samps = colnames(sorted_comb[, which(sorted_comb[2, ] == 'hubal_2016')])
gse83352_all_pdata = gse83352$GSE83352_series_matrix.txt.gz@phenoData@data

# Get row indices for relevant samples
indices = sapply(gse83352_samps, function(x) {which(gse83352_all_pdata$geo_accession == x)[1]})

gse83352_age = as.numeric(gse83352_all_pdata$`age:ch1`[indices])
names(gse83352_age) = rownames(gse83352_all_pdata[indices, ])

#### GSE87105 - no actual age! ####
# young: 31.2±0.86 y
# middle: 46.0±0.71 y
# old: 74.0±4.38 y

gse87105 = getGEO('GSE87105')
View(gse87105$GSE87105_series_matrix.txt.gz@phenoData@data)
gse87105_age = as.factor(gse87105$GSE87105_series_matrix.txt.gz@phenoData@data$`age:ch1`)
levels(gse87105_age) = c(46.0, 74.0, 31.2)
gse87105_age = as.numeric(as.character(gse87105_age))
names(gse87105_age) = rownames(gse87105$GSE87105_series_matrix.txt.gz@phenoData@data)

#### GSE98613 - no actual age! ####
# mean: 68.6 ± 16.8 y
# young: 24–42 (4) -- 33
# middle: 61–67 (5) -- 64
# old: 73–84 (15) -- 78.5

gse98613 = getGEO('GSE98613')
View(gse98613$GSE98613_series_matrix.txt.gz@phenoData@data)
gse98613_age = as.factor(gse98613$GSE98613_series_matrix.txt.gz@phenoData@data$`age:ch1`)
levels(gse98613_age) = c(64, 78.5, 33)
gse98613_age = as.numeric(as.character(gse98613_age))
names(gse98613_age) = rownames(gse98613$GSE98613_series_matrix.txt.gz@phenoData@data)

# Remove the repeat sample, 'GSM2601747'
gse98613_age = gse98613_age[!names(gse98613_age) == 'GSM2601747']

# Need to change to dataframe so that rbinding allows matching with names
all_ages = t(as.data.frame(c(gse40645_age, gse83352_age, gse87105_age, gse98613_age)))

#### Add ages to sorted comb table ####

sorted_comb_age = rbind(sorted_comb[1:4, ], all_ages, sorted_comb[6:nrow(sorted_comb), ])
rownames(sorted_comb_age)[5] = 'actual_age'

# Save new dataset

write.table(sorted_comb_age, file = 'sorted_comb_data_age.txt', sep = '\t', row.names = T, col.names = T)

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001")

sorted_comb_age = read.table('sorted_comb_data_age.txt', sep = '\t', stringsAsFactors = F)

#### Classify samples according to age ####

classify_age = function(dataset, young_thres, old_thres){
  
  young_count = 0; old_count = 0
  
  for (i in 1:ncol(dataset)){ # Loop through all columns
    
    if (dataset[5, i] <= young_thres) { # If as number and <= 50
      class_name = c('young')
      young_count = young_count + 1
    }
    else if (dataset[5, i] > old_thres) { # If as number and > 50
      class_name = c('old')
      old_count = old_count + 1
    }
    
    else {
      class_name = c('none')
    }
    
    dataset[5, i] = class_name # Substitute with either 'young' or 'old'
    
    i = i + 1
    
  }
  
  dataset = dataset[, !dataset[5, ] == 'none']
  
  #View(dataset)
  
  rownames(dataset)[5] = 'age'
  
  print(paste('old count:', old_count, sep = ' '))
  print(paste('young count:', young_count, sep = ' '))
  
  return(dataset)
  
}

# 1) Thres: <= 45, > 45

dataset_4545 = classify_age(sorted_comb_age, 45, 45)


# 2) Thres: <= 50, > 50

dataset_5050 = classify_age(sorted_comb_age, 50, 50)

write.table(dataset_5050, 'dataset_5050.txt', sep = '\t', row.names = T, col.names = T, quote = F)

# Only expression data
df_exp = sapply(dataset_5050[7:nrow(dataset_5050), ], as.numeric) # Get only expression data
rownames(df_exp) = rownames(dataset_5050[7:nrow(dataset_5050), ]) # Assign gene names as rownames

write.table(df_exp, 'ori.txt', sep = '\t', row.names = T, col.names = T, quote = F)

# 3) Thres: <= 55, > 55

dataset_5555 = classify_age(sorted_comb_age, 55, 55)
write.table(dataset_5555, 'dataset_5555.txt', sep = '\t', row.names = T, col.names = T, quote = F)

# 4) Thres: <= 40, > 40

dataset_6060 = classify_age(sorted_comb_age, 60, 60)
write.table(dataset_6060, 'dataset_6060.txt', sep = '\t', row.names = T, col.names = T, quote = F)
