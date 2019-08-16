################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Outlier Analysis


library(blackSheepR)

library(curl)
annotationtable = read.table(curl(paste0("https://raw.githubusercontent.com/",
                    "ruggleslab/blacksheep_supp/dev/vignettes/brca/",
                    "annotations_common_samples.csv")), header = TRUE,
                    row.names = 1, na.strings = c("", " ", "NA"), sep = ",",
                    check.names = FALSE, stringsAsFactors = FALSE)
colnames(annotationtable) = gsub(" ", "_", colnames(annotationtable))
compcols = annotationtable[,c("PAM50", "ER_Status", "PR_Status",
                        "GATA3_Mutation", "PIK3CA_Mutation", "TP53_Mutation")]

rnatable = read.table(curl(paste0("https://raw.githubusercontent.com/",
                    "ruggleslab/blacksheep_supp/dev/vignettes/brca/",
                    "rna_common_samples_data.csv")), header = TRUE,
                    row.names = 1, sep = ",", quote = "", check.names = FALSE)
detach("package:curl", unload=TRUE)


# Example Workflow
w - RNA
# In the following section - we will go through an example of using outlier
# analysis using RNA data. The inputted data is being supplied from
# [Github](https://github.com/ruggleslab/blacksheep_supp/tree/master) and is
# from breast cancer data from TCGA and CPTAC.

## Read in Annotation Data
# annotationfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
#                 "outlier-tool/data/brca/annotations_common_samples.csv")
# annotationtable = read.table(annotationfile, header = TRUE, row.names = 1,
#                 na.strings = c("", " ", "NA"), sep = ",", check.names = FALSE,
#                 stringsAsFactors = FALSE)
# colnames(annotationtable) = gsub(" ", "_", colnames(annotationtable))
# compcols = annotationtable[,c("PAM50", "ER_Status", "PR_Status",
#                         "GATA3_Mutation", "PIK3CA_Mutation", "TP53_Mutation")]

## FORMAT our annotation table
compcols[,4] = ifelse(is.na(compcols[,4]), "None", "Mutant")
compcols[,5] = ifelse(is.na(compcols[,5]), "None", "Mutant")
compcols[,6] = ifelse(is.na(compcols[,6]), "None", "Mutant")

expanded_compcols = make_comparison_columns(compcols[,1,drop=FALSE])

comptable = cbind(expanded_compcols, compcols[2:ncol(compcols)])

## Read in Count Data
rnacountfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                    "data/brca/rna_common_samples_data.csv")
rnatable = read.table(rnacountfile, header = TRUE, row.names = 1, sep = ",",
quote = "", check.names = FALSE)

outfilepath = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                "output/vignette_example/")
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

comptable[1:5,1:5]
dim(comptable)
rnatable[1:5,1:5]
dim(rnatable)

#### Create groupings
groupings = comparison_groupings(comptable)
## Print out the first 6 samples in each of our first 5 groupings
lapply(groupings, head)[1:5]

#### Make Outlier table
reftable_function_out = make_outlier_table(rnatable)
## See the names of the outputted objects
names(reftable_function_out)
## Assign them to individual variables
outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
sampmedtab = reftable_function_out$sampmedtab

## Note we will only use the outlier table - which looks like this now
outliertab[1:5,1:5]

#### Tabulate Outliers
grouptablist = count_outliers(groupings, outliertab)$grouptablist
lapply(grouptablist, head)[1:5]

#### Run Outlier Analysis
outlier_analysis_out = outlier_analysis(grouptablist)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1:5]

for (analysisnum in seq_len(length(outlier_analysis_out))) {
    subanalysis = outlier_analysis_out[[analysisnum]]
    analysisoutfile = paste0(outfilepath, "rna_",
                        names(outlier_analysis_out[analysisnum]), ".csv")
    write.table(subanalysis, analysisoutfile, row.names = FALSE,
                col.names = TRUE, sep = ",")
}

#### Plot Results using Heatmap Generating Function
hm_annotations = comptable

hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out,
        counttab = rnatable, metatable = hm_annotations, fdrcutoffvalue = 0.1)
hm1

hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out,
        counttab = rnatable, metatable = hm_annotations, fdrcutoffvalue = 0.1,
        write_out_plot = TRUE, outfilepath =
            paste0(outfilepath, "test_rna_hm1_"))

## To output heatmap to pdf outside of the function (not recommended)
# pdf(paste0(outfilepath, "test_rna_hm1.pdf"))
# hm1
# junk<-dev.off()




# Example Workflow - Phosphoprotein
# In the following section - we will go through an example of using outlier
# analysis using RNA data. The inputted data is being supplied from
# [Github](https://github.com/ruggleslab/blacksheep_supp/tree/master) and is
# from breast cancer data from TCGA and CPTAC.

#### Read in annotation Data
annotationfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                    "outlier-tool/data/brca/annotations_common_samples.csv")
annotationtable = read.table(annotationfile, header = TRUE, row.names = 1,
                na.strings = c("", " ", "NA"), sep = ",", check.names = FALSE)
compcols = annotationtable[,c("PAM50", "ER Status", "PR Status",
                        "GATA3 Mutation", "PIK3CA Mutation", "TP53 Mutation")]

## FORMAT our annotation table
compcols[,4] = ifelse(is.na(compcols[,4]), "None", "Mutant")
compcols[,5] = ifelse(is.na(compcols[,5]), "None", "Mutant")
compcols[,6] = ifelse(is.na(compcols[,6]), "None", "Mutant")

expanded_compcols = make_comparison_columns(compcols[,1,drop=FALSE])

comptable = cbind(expanded_compcols, compcols[2:ncol(compcols)])

## Read in Phosphoprotein Data
phosphocountfile = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                "outlier-tool/data/brca/phospho_common_samples_data.csv")
phosphotable = read.table(phosphocountfile, header = TRUE, row.names = 1,
                            sep = ",", quote = "", check.names = FALSE)

outfilepath = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                "output/vignette_example/")
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

comptable[1:5,1:5]
dim(comptable)
phosphotable[1:5,1:5]
dim(rnatable)
print("Testing for existance of outfilepath; <dir.exists(outfilepath)>")
dir.exists(outfilepath)


#### Create groupings
groupings = comparison_groupings(comptable)
## Print out the first 6 samples in each of our first 5 groupings
lapply(groupings, head)[1:5]

#### Make Outlier table
## Perform the function
# reftable_function_out = make_outlier_table(phosphotable,
#                         aggregate_features = TRUE, feature_delineator = "-")
reftable_function_out = make_outlier_table(phosphotable,
                                        analyze_negative_outliers = FALSE)
## See the names of the outputted objects
names(reftable_function_out)
## Assign them to individual variables
outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
sampmedtab = reftable_function_out$sampmedtab

## Note we will only use the outlier table - which looks like this now
outliertab[1:5,1:5]
aggposoutlierstab[1:5,1:5]
aggposfractiontab[1:5,1:5]

#### Tabulate Outliers
#grouptablist = count_outliers(groupings, aggposoutlierstab)
count_outliers_out = count_outliers(groupings, outliertab,
                        aggregate_features = TRUE, feature_delineator = "-")
grouptablist = count_outliers_out$grouptablist
aggoutliertab = count_outliers_out$aggoutliertab
aggfractiontab = count_outliers_out$aggfractiontab

lapply(grouptablist, head)[1:5]

#### Run Outlier Analysis
outlier_analysis_out = outlier_analysis(grouptablist,
                                        fraction_table = aggfractiontab)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1:5]

for (analysisnum in seq_len(length(outlier_analysis_out))) {
    subanalysis = outlier_analysis_out[[analysisnum]]
    analysisoutfile = paste0(outfilepath, "phospho_",
                            names(outlier_analysis_out[analysisnum]), ".csv")
    write.table(subanalysis, analysisoutfile, row.names = FALSE,
                col.names = TRUE, sep = ",")
}

#### Plot Results using Heatmap Generating Function
# hm_annotations = annotationtable[,c("PAM50", "ER Status", "PR Status",
#         "HER2 Status", "GATA3 Mutation", "PIK3CA Mutation", "TP53 Mutation")]
hm_annotations = comptable
hm_annotations = hm_annotations[order(hm_annotations[,1], hm_annotations[,2],
                    hm_annotations[,3] ,hm_annotations[,4], na.last = TRUE),]
hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out,
        counttab = aggfractiontab, metatable = hm_annotations,
        fdrcutoffvalue = 0.1)
hm1

hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out,
        counttab = aggfractiontab, metatable = hm_annotations,
        fdrcutoffvalue = 0.1, write_out_plot = TRUE, outfilepath =
                          paste0(outfilepath, "test_phopho_hm1_"))

## To output heatmap to pdf outside of the function (not recommended)
# pdf(paste0(outfilepath, "test_phospho_hm1.pdf"))
# hm1
# junk<-dev.off()


