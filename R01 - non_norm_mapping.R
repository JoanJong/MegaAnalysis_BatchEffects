## 2/7/19
## This R script standardises the format of non-normalised data from GSE98613, GSE87105, GSE83352 and GSE40645. 
## It also maps the probe IDs to the respective gene symbols.

library(GEOquery)
library(illuminaHumanv4.db)

setwd("C:/Users/User/OneDrive - Nanyang Technological University/Y2 S2/BS9001 Research Experience/R/win-library/aging_datasets")

## ------- IMPORT AND CLEAN NON-NORM DATA ------- ##

imported_98613 = read.table('GSE98613_non-normalized.txt', sep='\t', stringsAsFactors = FALSE)
imported_87105 = read.table('GSE87105_non-normalized.txt', sep='\t', stringsAsFactors = FALSE)
imported_83352 = read.table('GSE83352_non-normalized.txt', sep='\t', stringsAsFactors = FALSE)
imported_40645 = read.table('GSE40645_non-normalized.txt', sep='\t', stringsAsFactors = FALSE)

# Get the probe IDs
probe_98613 = imported_98613[4:nrow(imported_98613),1]
probe_87105 = imported_87105[2:nrow(imported_87105),1]
probe_83352 = imported_83352[2:nrow(imported_83352),1]
probe_40645 = imported_40645[2:nrow(imported_40645),1]

# Remove detection p-val columns - select only odd numbered columns
exp_98613 = imported_98613[3:nrow(imported_98613), c(FALSE, TRUE)]
exp_87105_1 = imported_87105[, c(FALSE, TRUE)]
exp_87105 = exp_87105_1[, 1:16] # Columns 18 onwards have entirely NA values and the experiment has only 16 samples anyway.

# GSE83352 has mismatched columns
tokeep_83352 = c(seq(2, 54, 2), 55, seq(58, 168, 2))
exp_83352 = subset(imported_83352, select = tokeep_83352) # All signal intensities, even if wrongly labelled
correct_names_idx_83352 = grep('E0', imported_83352[1, ])
correct_names_83352 = imported_83352[1, correct_names_idx_83352]
exp_83352[1, ] = correct_names_83352

exp_40645 = imported_40645[, c(FALSE, TRUE)]

# Get sample title
samp_98613 = exp_98613[1, ] # 24 samples
samp_87105 = exp_87105[1, ] # 16 samples
samp_83352 = exp_83352[1, ] # 84 samples
samp_40645 = exp_40645[1, ] # 29 samples

# Removing irrelevant columns and rows (sample titles)
nume_98613 = exp_98613[2:nrow(exp_98613),]
nume_87105 = exp_87105[2:nrow(exp_87105),]
nume_83352 = exp_83352[2:nrow(exp_83352),]
nume_40645 = exp_40645[2:nrow(exp_40645),]

# Making rownames the probe ID
rownames(nume_98613) = probe_98613
rownames(nume_87105) = probe_87105
rownames(nume_83352) = probe_83352
rownames(nume_40645) = probe_40645

# Making colnames the sample title
colnames(nume_98613) = samp_98613
colnames(nume_87105) = samp_87105
colnames(nume_83352) = samp_83352
colnames(nume_40645) = samp_40645

## -------- GETTING GENE SYMBOLS -------- ##

# Get a df of gene symbols matching to probe ID. DO NOT USE envir = illuminaHumanv4SYMBOL! (Many missing values)
gsym_98613 = data.frame(Gene = unlist(mget(x = as.character(probe_98613), envir = illuminaHumanv4SYMBOLREANNOTATED)))
gsym_87105 = data.frame(Gene = unlist(mget(x = as.character(probe_87105), envir = illuminaHumanv4SYMBOLREANNOTATED)))
gsym_83352 = data.frame(Gene = unlist(mget(x = as.character(probe_83352), ifnotfound = list(NA), envir = illuminaHumanv4SYMBOLREANNOTATED)))
gsym_40645 = data.frame(Gene = unlist(mget(x = as.character(probe_40645), ifnotfound = list(NA), envir = illuminaHumanv4SYMBOLREANNOTATED)))

# Can't get gene symbol for probe ILMN_1815977 (row 28758), so delete
# Use argument ifnotfound to skip the probes not in library!!

## ------- MAPPING GENE SYMBOLS TO PROBE IDS ------- ##

bind_98613 = cbind(gsym_98613, nume_98613) # 3576 missing values
bind_87105 = cbind(gsym_87105, nume_87105) # 3668 missing values
bind_83352 = cbind(gsym_83352, nume_83352) # 3034 missing values
bind_40645 = cbind(gsym_40645, nume_40645) # 11748 missing values

## ------- IMPUTATION ------- ##

# Remove rows that have NA in gene symbol column.

# Return only rows that don't have NA in gene symbol column
comp_98613 = bind_98613[!is.na(bind_98613$Gene), ] # 43655
comp_87105 = bind_87105[!is.na(bind_87105$Gene), ] # 43654
comp_83352 = bind_83352[!is.na(bind_83352$Gene), ] # 36236
comp_40645 = bind_40645[!is.na(bind_40645$Gene), ] # 37055

# Check if there are missing values in the expression data columns.
any(is.na(comp_98613)) # FALSE
any(is.na(comp_87105)) # FALSE
any(is.na(comp_83352))
any(is.na(comp_40645))

## ------- EXPORT ------- ##

write.table(comp_98613, file='./nn_mapped_txt/GSE98613_nn_mapped.txt', sep="\t", row.names=T, col.names=T, quote=F)
write.table(comp_87105, file='./nn_mapped_txt/GSE87105_nn_mapped.txt', sep="\t", row.names=T, col.names=T, quote=F)
write.table(comp_83352, file='./nn_mapped_txt/GSE83352_nn_mapped.txt', sep="\t", row.names=T, col.names=T, quote=F)
write.table(comp_40645, file='./nn_mapped_txt/GSE40645_nn_mapped.txt', sep="\t", row.names=T, col.names=T, quote=F)

