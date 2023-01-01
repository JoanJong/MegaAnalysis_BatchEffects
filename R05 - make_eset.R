## 8/7/20

## Make eSet using phenodata from dataset_5050.txt and non-normliased data from ori.txt

#### --- 1) Import exprs and pheno data --- ####

## Import exprs data (ori)

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001/normalised_datasets_50")
exprs_imp = read.table('ori.txt', sep = '\t', header = T, stringsAsFactors = F)

## Log2-transform exprs data

exprs_raw = as.matrix(exprs_imp)
exprs_imp = as.matrix(log2(exprs_imp))

## Import pheno data

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y3 S1/CY2001 Research Attachment/Codes and Data/CY2001/different_age_thresholds/50")
pheno_imp = read.table('dataset_5050.txt', sep = '\t', header = T, stringsAsFactors = FALSE)
pheno_imp = as.data.frame(t(pheno_imp[1:6, ]))

levels(pheno_imp$submission_date) = c('gonzalez', 'hubal', 'mercken', 'raz')

colnames(pheno_imp) = c('geo_accession', 'batch_id', 'platform_id', 'sex', 'age', 'tissue_type')

#### --- 2) Checking if exprs and pheno data are compatible --- ####

## Check if order of samples is the same
# Exprs: rownames are genes, colnames are samples
# Pheno: rownames are samples, colnames are pheno variables

all(colnames(exprs_imp) == rownames(pheno_imp)) # TRUE

#### --- 3) Making the eSet --- ####

## Need to make sure phenodata is an annotated dataframe

phenodata = new('AnnotatedDataFrame',
                data = pheno_imp)

new_eset = ExpressionSet(assayData = exprs_imp,
                         phenoData = phenodata)

#### --- 4) Exporting eSet --- ####

setwd("C:/Users/User/OneDrive - Nanyang Technological University(1)/Y4 S1 (OFYP)/Codes/exprs_data")
saveRDS(file = 'ori_50_eset.rds', new_eset)