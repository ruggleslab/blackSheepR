
## Read in the data
annotationfile = "/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/data/brca/annotations_common_samples.csv"
annotationtable = read.table(annotationfile, header = TRUE, row.names = 1, na.strings = c("", " ", "NA"), sep = ",")
comptable = annotationtable[,(ncol(annotationtable)-4):ncol(annotationtable)]

rnacountfile = "/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/data/brca/rna_common_samples_data.csv"
rnatable = read.table(rnacountfile, header = TRUE, row.names = 1, sep = ",", quote = "", check.names = FALSE)

outfilepath = "/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/output_VIGNETTE/"
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

comptable[1:5,1:5]
dim(comptable)
rnatable[1:5,1:5]
dim(rnatable)
print("Testing for existance of outfilepath; <dir.exists(outfilepath)>")
dir.exists(outfilepath)


## Create groupings
groupings = comparison_groupings(comptable)
## Print out the first 6 samples in each of our first 5 groupings
lapply(groupings, head)[1:5]


## Make the outlier table
reftable_function_out = make_outlier_table(rnatable)
## See the names of the outputted objects
names(reftable_function_out)
## Assign them to individual variables
outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
sampmedtab = reftable_function_out$sampmedtab

## Note we will only use the outlier table - which looks like this now
outliertab[1:5,1:5]


## Tabulate number of outliers
grouptablist = count_outliers(groupings, outliertab)
lapply(grouptablist, head)[1:5]


## Perfrom the outlier abalysis
outlier_analysis_out = outlier_analysis(grouptablist)
names(outlier_analysis_out)
lapply(outlier_analysis_out, head)[1:5]
