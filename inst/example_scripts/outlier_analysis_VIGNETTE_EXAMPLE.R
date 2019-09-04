################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Outlier Analysis

## Installation:
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("blacksheepr")
library(blacksheepr)


## Input Data
library(blacksheepr)
data("sample_values")
sample_values[1:5,1:5]

##### 2 - Annotation data
data("sample_annotations")
sample_annotations[1:5,1:5]


# Example Workflow - RNA

### Read in Annotation
library(curl)
annotationtable = read.table(curl(paste0("https://raw.githubusercontent.com/",
                        "ruggleslab/blacksheep_supp/dev/vignettes/brca/",
                        "annotations_common_samples.csv")), header = TRUE,
                        row.names = 1, na.strings = c("", " ", "NA"), sep = ",",
                        check.names = FALSE, stringsAsFactors = FALSE)
colnames(annotationtable) = gsub(" ", "_", colnames(annotationtable))
compcols = annotationtable[,c("PAM50", "ER_Status", "PR_Status",
                        "GATA3_Mutation", "PIK3CA_Mutation", "TP53_Mutation")]

outfilepath = getwd()

compcols[1:5,1:5]
dim(compcols)

#### Formatting our annotation data
compcols[,"GATA3_Mutation"] = ifelse(is.na(compcols[,4]), "None", "Mutant")
compcols[,"PIK3CA_Mutation"] = ifelse(is.na(compcols[,5]), "None", "Mutant")
compcols[,"TP53_Mutation"] = ifelse(is.na(compcols[,6]), "None", "Mutant")
head(compcols)

## Use the make_comparison_columns function to create binary columns
expanded_compcols = make_comparison_columns(compcols[,1,drop=FALSE])
## Append new columns to the comparison annotation table
comptable = cbind.data.frame(expanded_compcols, compcols[2:ncol(compcols)],
                            stringsAsFactors = FALSE)

### Reading in the RNA data
rnatable = read.table(curl(paste0("https://raw.githubusercontent.com/",
                    "ruggleslab/blacksheep_supp/dev/vignettes/brca/",
                    "rna_common_samples_data.csv")), header = TRUE,
                    row.names = 1, sep = ",", quote = "", check.names = FALSE)

rnatable[1:5,1:5]
dim(rnatable)


### Creating a SummarizedExperiment from our data.
suppressPackageStartupMessages(library(SummarizedExperiment))

blacksheep_SE = SummarizedExperiment(
    assays=list(counts=as.matrix(rnatable)),
    colData=DataFrame(comptable))
blacksheep_SE

### Running blacksheep
deva_out = deva(blacksheep_SE,
                analyze_negative_outliers = FALSE, aggregate_features = FALSE,
                fdrcutoffvalue = 0.1, write_out = FALSE)

names(deva_out$pos_outlier_analysis)
head(deva_out$pos_outlier_analysis$
        outlieranalysis_for_PAM50_Her2__Her2_vs_PAM50_Her2__not_Her2)
deva_out$significant_pos_heatmaps$
    print_outlieranalysis_for_PAM50_Basal__not_Basal_vs_PAM50_Basal__Basal

## NOT RUN
## To output separately to pdf
#pdf(outfile.pdf)
#deva_out$significant_pos_heatmaps$
#    print_outlieranalysis_for_PAM50_Basal__not_Basal_vs_PAM50_Basal__Basal
#dev.off()


# Piecewise analysis - RNA
groupings = comparison_groupings(comptable)
## Print out the first 6 samples in each of our first 5 groupings
lapply(groupings, head)[1:5]


### Make Outlier table
reftable_function_out = make_outlier_table(rnatable,
                                analyze_negative_outliers = FALSE)
## See the names of the outputted objects
names(reftable_function_out)
## Assign them to individual variables
outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
sampmedtab = reftable_function_out$sampmedtab

## Note we will only use the outlier table - which looks like this now
outliertab[1:5,1:5]


### Tabulate Outliers
count_outliers_out = count_outliers(groupings, outliertab)
grouptablist = count_outliers_out$grouptablist

names(grouptablist)
names(grouptablist$PAM50_Her2__Her2)
head(grouptablist$PAM50_Her2__Her2$feature_counts)
grouptablist$PAM50_Her2__Her2$samples



### Run Outlier Analysis
outlier_analysis_out = outlier_analysis(grouptablist)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1]

### Plot Results using Heatmap Generating Function
## This will sort every column in <comptable> sequentially before running
plottable = comptable[do.call(order, c(decreasing = TRUE,
data.frame(comptable[,1:ncol(comptable)]))),]
hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out,
                    counttab = rnatable, metatable = plottable,
                    fdrcutoffvalue = 0.1)
hm1[[4]]

## To output heatmap to pdf outside of the function
#pdf(paste0(outfilepath, "test_hm1.pdf"))
#hm1
#junk<-dev.off()






# Example Workflow - Phosphoprotein
# In the following section - we will go through an example of using outlier
# analysis using phospho data. The inputted data is being supplied from
# [Github](https://github.com/ruggleslab/blacksheep_supp/tree/master) and is
# from breast cancer data from TCGA and CPTAC.

#### Read in annotation Data
library(curl)
annotationtable = read.table(curl(paste0("https://raw.githubusercontent.com/",
                        "ruggleslab/blacksheep_supp/dev/vignettes/brca/",
                        "annotations_common_samples.csv")), header = TRUE,
                        row.names = 1, na.strings = c("", " ", "NA"), sep = ",",
                        check.names = FALSE, stringsAsFactors = FALSE)
colnames(annotationtable) = gsub(" ", "_", colnames(annotationtable))
compcols = annotationtable[,c("PAM50", "ER_Status", "PR_Status",
                        "GATA3_Mutation", "PIK3CA_Mutation", "TP53_Mutation")]
## FORMAT our annotation table
compcols[,4] = ifelse(is.na(compcols[,4]), "None", "Mutant")
compcols[,5] = ifelse(is.na(compcols[,5]), "None", "Mutant")
compcols[,6] = ifelse(is.na(compcols[,6]), "None", "Mutant")

expanded_compcols = make_comparison_columns(compcols[,1,drop=FALSE])

comptable = cbind(expanded_compcols, compcols[2:ncol(compcols)])

## Read in Phosphoprotein Data
phosphotable = read.table(curl(paste0("https://raw.githubusercontent.com/",
                    "ruggleslab/blacksheep_supp/dev/vignettes/brca/",
                    "phospho_common_samples_data.csv")), header = TRUE,
                    row.names = 1, sep = ",", quote = "", check.names = FALSE)

outfilepath = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
                "output/vignette_example/phosphos_noagg_neg/")
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
reftable_function_out = make_outlier_table(intable = phosphotable,
                                        analyze_negative_outliers = TRUE)
## See the names of the outputted objects
names(reftable_function_out)
## Assign them to individual variables
outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
sampmedtab = reftable_function_out$sampmedtab
lowerboundtab = reftable_function_out$lowerboundtab

## Note we will only use the outlier table - which looks like this now
outliertab[1:5,1:5]

#### Tabulate Outliers
#grouptablist = count_outliers(groupings, aggposoutlierstab)
count_outliers_out = count_outliers(groupings = groupings,
                        outliertab = outliertab,
                        aggregate_features = FALSE, feature_delineator = "-")
grouptablist = count_outliers_out$grouptablist
aggoutliertab = count_outliers_out$aggoutliertab
fractiontab = count_outliers_out$fractiontab

lapply(grouptablist, head)[1:5]

#### Run Outlier Analysis
outlier_analysis_out = outlier_analysis(
    grouptablist = grouptablist,
    fraction_table = fractiontab, fraction_samples_cutoff = 0.3,
    write_out_tables = TRUE, outfilepath = outfilepath)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1:5]

for (analysisnum in seq_len(length(outlier_analysis_out))) {
    subanalysis = outlier_analysis_out[[analysisnum]]
    analysisoutfile = paste0(outfilepath, "manual_phospho_",
                            names(outlier_analysis_out[analysisnum]), ".csv")
    write.table(subanalysis, analysisoutfile, row.names = FALSE,
                col.names = TRUE, sep = ",")
}

#### Plot Results using Heatmap Generating Function
# hm_annotations = annotationtable[,c("PAM50", "ER Status", "PR Status",
#         "HER2 Status", "GATA3 Mutation", "PIK3CA Mutation", "TP53 Mutation")]
plottable = comptable[do.call(order, c(decreasing = TRUE,
                                data.frame(comptable[,1:ncol(comptable)]))),]
hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out,
        counttab = fractiontab, metatable = plottable,
        fdrcutoffvalue = 0.01)
hm1

hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out,
        counttab = fractiontab, metatable = plottable,
        fdrcutoffvalue = 0.01, write_out_plot = TRUE, outfilepath =
                          paste0(outfilepath, "test_phopho_hm1_"))

## To output heatmap to pdf outside of the function (not recommended)
pdf(paste0(outfilepath, "manual_test_phospho_hm1.pdf"))
hm1
junk<-dev.off()


# ### Running blacksheep
# The blacksheep function has a number of steps to it that are individually
# described below. Note though that the individual steps only need to be used
# for specific query or alteration. In the general case, the `deva` function on
# its own should be sufficient for the desired analysis.
suppressPackageStartupMessages(library(SummarizedExperiment))

deva_SE = SummarizedExperiment(
    assays=list(counts=as.matrix(phosphotable)),
    colData=DataFrame(comptable))


dir.create(paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/",
           "output/vignette_example/deva_phosphos_noagg_neg/"),
           recursive = TRUE, showWarnings = FALSE)
deva_out = deva(se = deva_SE,
        analyze_negative_outliers = TRUE, aggregate_features = FALSE,
        feature_delineator = "-", fraction_samples_cutoff = 0.3,
        fdrcutoffvalue = 0.01, write_out = TRUE,
        outfilepath = paste0("/Users/tosh/Desktop/Ruggles_Lab/projects/",
                             "outlier-tool/output/vignette_example/",
                             "deva_phosphos_noagg_neg/"))


