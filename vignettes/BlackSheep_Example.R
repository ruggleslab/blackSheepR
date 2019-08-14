## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install package, eval=FALSE-----------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("blackSheepR")

## ----library call--------------------------------------------------------
library(blackSheepR)

## ----countdata example, echo=FALSE---------------------------------------
library(blackSheepR)
data("sample_values")
sample_values[1:5,1:5]

## ----annotation example, echo=FALSE--------------------------------------
data("sample_annotations")
sample_annotations[1:5,1:5]

## ----read in data - rna, echo = FALSE------------------------------------
annotationfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                        "outlier-tool/data/brca/annotations_common_samples.csv")
annotationtable = read.table(annotationfile, header = TRUE, row.names = 1, 
                na.strings = c("", " ", "NA"), sep = ",", check.names = FALSE)
comptable = annotationtable[,(ncol(annotationtable)-4):ncol(annotationtable)]
colnames(comptable) = c("PIK3CA_helical_mutant", "PIK3CA_kinase_mutant", 
                        "TP53_Nonsense", "TP53_Missense_all", 
                        "TP53_Missense_DNABD")

rnacountfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                      "data/brca/rna_common_samples_data.csv")
rnatable = read.table(rnacountfile, header = TRUE, row.names = 1, sep = ",", 
                      quote = "", check.names = FALSE)

outfilepath = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                     "output_VIGNETTE/")
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

comptable[1:5,1:5]
dim(comptable)
rnatable[1:5,1:5]
dim(rnatable)
print("Testing for existance of outfilepath; <dir.exists(outfilepath)>")
dir.exists(outfilepath)


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
grouptablist = count_outliers(groupings, outliertab)
lapply(grouptablist, head)[1:5]

## ----outlier analysis - rna----------------------------------------------
outlier_analysis_out = outlier_analysis(grouptablist)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1:5]

## ----heatmap plotting - rna----------------------------------------------
hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out, 
                      counttab = rnatable, metatable = annotationtable, 
                      fdrcutoffvalue = 0.1)
hm1

## To output heatmap to pdf outside of the function
#pdf(paste0(outfilepath, "test_hm1.pdf"))
#hm1
#junk<-dev.off()

## ----read in data - phospho, echo = FALSE--------------------------------
annotationfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                        "outlier-tool/data/brca/annotations_common_samples.csv")
annotationtable = read.table(annotationfile, header = TRUE, row.names = 1, 
                na.strings = c("", " ", "NA"), sep = ",", check.names = FALSE)
comptable = annotationtable[,(ncol(annotationtable)-4):ncol(annotationtable)]
colnames(comptable) = c("PIK3CA_helical_mutant", "PIK3CA_kinase_mutant", 
                        "TP53_Nonsense", "TP53_Missense_all", 
                        "TP53_Missense_DNABD")

phosphocountfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                    "outlier-tool/data/brca/phospho_common_samples_data.csv")
phosphotable = read.table(phosphocountfile, header = TRUE, row.names = 1, 
                        sep = ",", quote = "", check.names = FALSE)

outfilepath = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                     "output_VIGNETTE/")
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

comptable[1:5,1:5]
dim(comptable)
phosphotable[1:5,1:5]
dim(rnatable)
print("Testing for existance of outfilepath; <dir.exists(outfilepath)>")
dir.exists(outfilepath)


## ----groupings - phospho-------------------------------------------------
groupings = comparison_groupings(comptable)
## Print out the first 6 samples in each of our first 5 groupings
lapply(groupings, head)[1:5]


## ----make outlier table - phospho----------------------------------------
## Perform the function
reftable_function_out = make_outlier_table(phosphotable, 
                            aggregate_features = TRUE, feature_delineator = "-")
## See the names of the outputted objects
names(reftable_function_out)
## Assign them to individual variables
outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
sampmedtab = reftable_function_out$sampmedtab
aggposoutlierstab = reftable_function_out$aggposoutlierstab
aggposfractiontab = reftable_function_out$aggposfractiontab

## Note we will only use the outlier table - which looks like this now
outliertab[1:5,1:5]
aggposoutlierstab[1:5,1:5]
aggposfractiontab[1:5,1:5]


## ----groupingtablist - phospho-------------------------------------------
grouptablist = count_outliers(groupings, aggposoutlierstab)
lapply(grouptablist, head)[1:5]

## ----outlier analysis - phospho------------------------------------------
outlier_analysis_out = outlier_analysis(grouptablist)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1:5]

## ----heatmap plotting - phospho------------------------------------------
hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out, 
                counttab = aggposfractiontab, metatable = annotationtable, 
                fdrcutoffvalue = 0.1)
hm1

## To output heatmap to pdf outside of the function
#pdf(paste0(outfilepath, "test_hm1.pdf"))
#hm1
#junk<-dev.off()

