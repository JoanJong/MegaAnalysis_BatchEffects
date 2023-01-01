
## 21/6/19
## This R script adds additional data to dataframes compressed by df_compress.R.
## The resultant dataframes are then combined together into one dataframe.
## The output from this script is used for batch effect analysis.
## Additional data to be added includes:
## - Study name
## - Platform (Affymetrix/Illumina)
## - Cell type
## - Class label
## - Sex
## - Age
## - Ethnicity
## - Intervention
## - Gene symbols

library(GEOquery)
library(dplyr)

# Set working directory

## ------- 1) IMPORT ------- ##

setwd("C:/Users/User/OneDrive - Nanyang Technological University/Y2 S2/BS9001 Research Experience/R/win-library/aging_datasets/nn_compressed_df")

# Dataframe of gene symbol against sample name, showing expression data. All gene symbols appear only once.
cps_98613 = read.table("GSE98613_nn_mapped_cps.txt", header=T, sep="\t", stringsAsFactors = FALSE) # 24 samples
cps_87105 = read.table("GSE87105_nn_mapped_cps.txt", header=T, sep="\t", stringsAsFactors = FALSE) # 16 samples
cps_83352 = read.table("GSE83352_nn_mapped_cps.txt", header=T, sep="\t", stringsAsFactors = FALSE) # 84 samples
cps_40645 = read.table("GSE40645_nn_mapped_cps.txt", header=T, sep="\t", stringsAsFactors = FALSE) # 29 samples

## ------- 2) GETTING GSM NAME ------- ##

# In the non-normalised txt file, no GSM labels are given. Need to retrieve this information using pData.

setwd("C:/Users/User/OneDrive - Nanyang Technological University/Y2 S2/BS9001 Research Experience/R/win-library/aging_datasets")

geo_98613 = getGEO(filename = 'GSE98613_series_matrix.txt.gz')
geo_87105 = getGEO(filename = 'GSE87105_series_matrix.txt.gz')
geo_83352 = getGEO(filename = 'GSE83352_series_matrix.txt.gz')
geo_40645 = getGEO(filename = 'GSE40645_series_matrix.txt.gz')

pheno_98613 = pData(geo_98613)
pheno_87105 = pData(geo_87105)
pheno_83352 = pData(geo_83352)
pheno_40645 = pData(geo_40645)

# Names of the columns to extract
colex_98613 = c('geo_accession', 'submission_date', 'platform_id', 'gender:ch1', 'age:ch1', 'source_name_ch1') 
colex_87105 = c('geo_accession', 'submission_date', 'platform_id', 'gender:ch1', 'age:ch1', 'tissue:ch1')
colex_83352 = c('geo_accession', 'submission_date', 'platform_id', 'gender:ch1', 'age:ch1', 'tissue:ch1')
colex_40645 = c('geo_accession', 'submission_date', 'platform_id', 'gender:ch1', 'age:ch1', 'source_name_ch1')

# Subset out relevant portion of pheno dataframe
phedf_98613 = t(subset(pheno_98613, select = colex_98613)) # matching
phedf_87105 = t(subset(pheno_87105, select = colex_87105)) # matching
phedf_83352 = t(subset(pheno_83352, select = colex_83352)) # not matching
phedf_40645 = t(subset(pheno_40645, select = colex_40645)) # matching

# Need to rearrange the samples in GSE83352!
vec_order = as.character(pheno_83352[, 1])
cps_83352 = cps_83352[, match(vec_order, colnames(cps_83352))] # match columns in cps_83352 to the order of sample titles in pheno_83352

# Set names of rows and columns
new_phenames = c('geo_accession', 'submission_date', 'platform_id', 'sex', 'age', 'tissue_type')

rownames(phedf_98613) = new_phenames
rownames(phedf_87105) = new_phenames
rownames(phedf_83352) = new_phenames
rownames(phedf_40645) = new_phenames

# Now all rownames are gene symbols and all col names are GSM names.
colnames(cps_98613) = phedf_98613[1, ]
colnames(cps_87105) = phedf_87105[1, ]
colnames(cps_83352) = phedf_83352[1, ]
colnames(cps_40645) = phedf_40645[1, ]

# Find intersecting genes
inter_40645_83352 = intersect(rownames(cps_40645), rownames(cps_83352))
inter_87105_98613 = intersect(rownames(cps_87105), rownames(cps_98613))
inter_all = intersect(inter_40645_83352, inter_87105_98613) # 21053 genes in common
num_98613 = cps_98613[inter_all,] 
num_87105 = cps_87105[inter_all,] 
num_83352 = cps_83352[inter_all,] 
num_40645 = cps_40645[inter_all,] 

# Bind expression data from 3 studies
comb_numeric_1 = cbind(num_40645, num_83352)
comb_numeric_2 = cbind(num_87105, num_98613)
comb_numeric_all = cbind(comb_numeric_1, comb_numeric_2) # 153 samples
comb_num_fac = t(apply(comb_numeric_all, 1, as.factor))

# Export numeric only table (unfiltered)
setwd("C:/Users/User/OneDrive - Nanyang Technological University/Y2 S2/BS9001 Research Experience/R/win-library/aging_datasets")
write.table(comb_numeric_all, file='./comb_num_unfilt.txt', sep="\t", row.names=T, col.names=T, quote=F)

# Bind pheno data to expression df
all_pheno_data_1 = cbind(phedf_40645, phedf_83352)
all_pheno_data_2 = cbind(phedf_87105, phedf_98613)
all_pheno_data = cbind(all_pheno_data_1, all_pheno_data_2)
all_data = as.data.frame(rbind(all_pheno_data, comb_num_fac), stringsAsFactors = FALSE) # GET RID OF FACTORS!!!

## ------- CLEANING ------- ##

# Some samples need to be removed.

# Delete repeat sample in GSE98613 - GSM2601747
# Delete samples that underwent interventions
# GSE87105 and GSE40645 are clean
clean_comb_98613 = subset(all_data, select = -GSM2601747) # Remove 1 column
gsmof83352 = phedf_83352[1,]
drop_83352 = gsmof83352[c(TRUE, FALSE)]
clean_comb_83352 = clean_comb_98613[, -which(colnames(clean_comb_98613) %in% drop_83352)] # Remove 42 columns
# Total samples available: 110
clean_comb = clean_comb_83352

# Classification of old, middle-aged and young different for 4 studies -- need to standardise
# > levels(as.factor(all_data[5,]))
# [1] "28"           "29"           "34"           "35"           "35 years old" "36 years old" "37 years old" "39 years old"
# [9] "40 years old" "41"           "41 years old" "42 years old" "43 years old" "44"           "44 years old" "45"          
# [17] "45 years old" "46"           "47"           "47 years old" "48"           "48 years old" "50"           "51"          
# [25] "52"           "53"           "54"           "55"           "56"           "56 years old" "57 years old" "58 years old"
# [33] "59"           "60"           "60 years old" "63"           "64"           "65"           "66 years old" "68 years old"
# [41] "73 years old" "77 years old" "85 years old" "87 years old" "89 years old" "Middle-aged"  "Old"          "Young"    

# Batch factors as submission date

i = 1
young_count = 0
old_count = 0

# Standardising class labels
for (class_name in clean_comb[5, ]) {
  age_char = substr(class_name, start = 1, stop = 2)
  age_num = as.numeric(age_char)
  print(age_num)
  
  # Check if age given as number or as string
  if (is.na(age_num)) { # If as string
    if (age_char == 'Yo'|age_char == 'Mi'){
      class_name = c('young')
      young_count = young_count + 1
    }
    else if (age_char == 'Ol'){
      class_name = c('old')
      old_count = old_count + 1
    }
  }
  
  else if (age_num <= 50) { # If as number and <= 50
    class_name = c('young')
    young_count = young_count + 1
  }
  else if (age_num > 50) { # If as number and > 50
    class_name = c('old')
    old_count = old_count + 1
  }
    
  clean_comb[5,i] = class_name # Substitute with either 'young' or 'old'
  
  i = i + 1
}

print(old_count)
print(young_count)

# Standardise study names

x = 1
study_authors = c('raz_2012', 'hubal_2016', 'mercken_2016',
                  'gonzalez_2017')
for (date in clean_comb[2, ]) {
  if (date == 'Sep 06 2012'){
    clean_comb[2, x] = study_authors[1]
  }
  else if (date == 'Jun 14 2016'){
    clean_comb[2, x] = study_authors[2]
  }
  else if (date == 'Sep 19 2016'){
    clean_comb[2, x] = study_authors[3]
  }
  else if (date == 'May 05 2017'){
    clean_comb[2, x] = study_authors[4]
  }
  x = x + 1
}

# Standardise gender
clean_comb[4, ] = tolower(clean_comb[4, ])

# Stadardise tissue -- all from vastus lateralis
clean_comb[6, ] = c('vastus_lateralis')

setwd("C:/Users/User/OneDrive - Nanyang Technological University/Y2 S2/BS9001 Research Experience/R/win-library/aging_datasets")
write.table(clean_comb, file='./sorted_comb_data.txt', sep="\t", row.names=T, col.names=T, quote=F)

