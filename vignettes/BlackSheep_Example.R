## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install package, eval=FALSE-----------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("blackSheepR")

## ----library call--------------------------------------------------------
library(blackSheepR)

## ----countdata example---------------------------------------------------
library(blackSheepR)
data("sample_values")
sample_values[1:5,1:5]

## ----annotation example--------------------------------------------------
data("sample_annotations")
sample_annotations[1:5,1:5]

## ----read in data - rna, echo = FALSE------------------------------------
annotationfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                        "outlier-tool/data/brca/annotations_common_samples.csv")
annotationtable = read.table(annotationfile, header = TRUE, row.names = 1, 
                na.strings = c("", " ", "NA"), sep = ",", check.names = FALSE)
colnames(annotationtable) = gsub(" ", "_", colnames(annotationtable))
compcols = annotationtable[,c("PAM50", "ER_Status", "PR_Status",
                        "GATA3_Mutation", "PIK3CA_Mutation", "TP53_Mutation")]

rnacountfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                      "data/brca/rna_common_samples_data.csv")
rnatable = read.table(rnacountfile, header = TRUE, row.names = 1, sep = ",", 
                      quote = "", check.names = FALSE)

outfilepath = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                     "output_VIGNETTE/")
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

compcols[1:5,1:5]
dim(compcols)
rnatable[1:5,1:5]
dim(rnatable)
print("Testing for existance of outfilepath; <dir.exists(outfilepath)>")
dir.exists(outfilepath)


## ----format annotation data1 - rna---------------------------------------
## FORMAT our annotation table
compcols[,4] = ifelse(is.na(compcols[,4]), "None", "Mutant")
compcols[,5] = ifelse(is.na(compcols[,5]), "None", "Mutant")
compcols[,6] = ifelse(is.na(compcols[,6]), "None", "Mutant")
head(compcols)

## ----format annotation data2 - rna---------------------------------------
## Use the make_comparison_columns function to create binary columns
expanded_compcols = make_comparison_columns(compcols[,1,drop=FALSE])
## Append new columns to the comparison annotation table
comptable = cbind(expanded_compcols, compcols[2:ncol(compcols)])
head(comptable)

## ----groupings - rna-----------------------------------------------------
groupings = comparison_groupings(comptable)
## Print out the first 6 samples in each of our first 5 groupings
lapply(groupings, head)[1:5]


## ----make outlier table - rna--------------------------------------------
## Perform the function
reftable_function_out = make_outlier_table(rnatable)
## See the names of the outputted objects
names(reftable_function_out)
## Assign them to individual variables
outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
sampmedtab = reftable_function_out$sampmedtab

## Note we will only use the outlier table - which looks like this now
outliertab[1:5,1:5]


## ----groupingtablist - rna-----------------------------------------------
count_outliers_out = count_outliers(groupings, outliertab)
grouptablist = count_outliers_out$grouptablist

names(grouptablist)

## ------------------------------------------------------------------------
names(grouptablist[[1]])

## ------------------------------------------------------------------------
head(grouptablist[[1]][[1]])

## ------------------------------------------------------------------------
grouptablist[[1]][[2]]

## ----outlier analysis - rna----------------------------------------------
## Run the outlier analysis function
outlier_analysis_out = outlier_analysis(grouptablist)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1]

## ----heatmap plotting - rna, fig.keep="none"-----------------------------
hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out, 
                      counttab = rnatable, metatable = comptable, 
                      fdrcutoffvalue = 0.1)

## To output heatmap to pdf outside of the function
#pdf(paste0(outfilepath, "test_hm1.pdf"))
#hm1
#junk<-dev.off()

## ----heatmap plotting 1example - rna-------------------------------------
hm1[[1]]

## ----read in data - phospho, echo = FALSE--------------------------------
annotationfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                        "outlier-tool/data/brca/annotations_common_samples.csv")
annotationtable = read.table(annotationfile, header = TRUE, row.names = 1, 
                na.strings = c("", " ", "NA"), sep = ",", check.names = FALSE)
colnames(annotationtable) = gsub(" ", "_", colnames(annotationtable))
compcols = annotationtable[,c("PAM50", "ER_Status", "PR_Status",
                        "GATA3_Mutation", "PIK3CA_Mutation", "TP53_Mutation")]

phosphocountfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                    "outlier-tool/data/brca/phospho_common_samples_data.csv")
phosphotable = read.table(phosphocountfile, header = TRUE, row.names = 1, 
                        sep = ",", quote = "", check.names = FALSE)

outfilepath = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                     "output_VIGNETTE/")
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

compcols[1:5,1:5]
dim(compcols)
rnatable[1:5,1:5]
dim(rnatable)
print("Testing for existance of outfilepath; <dir.exists(outfilepath)>")
dir.exists(outfilepath)


## ----format annotation data1 - phospho-----------------------------------
## FORMAT our annotation table
compcols[,4] = ifelse(is.na(compcols[,4]), "None", "Mutant")
compcols[,5] = ifelse(is.na(compcols[,5]), "None", "Mutant")
compcols[,6] = ifelse(is.na(compcols[,6]), "None", "Mutant")
head(compcols)

## ----format annotation data2 - phospho-----------------------------------
## Use the make_comparison_columns function to create binary columns
expanded_compcols = make_comparison_columns(compcols[,1,drop=FALSE])
## Append new columns to the comparison annotation table
comptable = cbind(expanded_compcols, compcols[2:ncol(compcols)])
head(comptable)

## ----groupings - phospho-------------------------------------------------
groupings = comparison_groupings(comptable)
## Print out the first 6 samples in each of our first 5 groupings
lapply(groupings, head)[1:5]

## ----make outlier table - phospho----------------------------------------
## Perform the function
reftable_function_out = make_outlier_table(phosphotable)
## See the names of the outputted objects
names(reftable_function_out)
## Assign them to individual variables
outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
sampmedtab = reftable_function_out$sampmedtab

## Note we will only use the outlier table - which looks like this now
outliertab[1:5,1:5]

## ----groupingtablist - phospho-------------------------------------------
count_outliers_out = count_outliers(groupings, outliertab, 
                        aggregate_features = TRUE, feature_delineator = "-")
grouptablist = count_outliers_out$grouptablist
aggoutliertab = count_outliers_out$aggoutliertab
aggfractiontab = count_outliers_out$aggfractiontab

names(grouptablist)

## ------------------------------------------------------------------------
names(grouptablist[[1]])

## ------------------------------------------------------------------------
head(grouptablist[[1]][[1]])

## ------------------------------------------------------------------------
grouptablist[[1]][[2]]

## ----outlier analysis - phospho------------------------------------------
outlier_analysis_out = outlier_analysis(grouptablist,
                                        fraction_table = aggfractiontab,
                                        fraction_samples_cutoff = 0.3)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1]

