################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Outlier Analysis

## Load in Libraries
print("Beginning Outlier Analysis")
packagelist = c("tools")
junk <- lapply(packagelist, function(xxx) suppressMessages(require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

#####

#normcounttabfile = "/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/data/new/sample_endo.csv"
#metatablefile = "/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/data/new/sample_annotations.csv"
#outfilepath = "/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/output/"

# sample_endo = t(read.table(normcounttabfile, header = TRUE, row.names = 1, sep = ",", check.names = FALSE))
# sample_annotations = read.table(metatablefile, header = TRUE, row.names = 1, sep = ifelse(tools::file_ext(metatablefile)=="txt", "\t", ","), stringsAsFactors = FALSE, check.names = FALSE, na.strings = c(NA, "NA", ""))
# namecheck = intersect(colnames(normcounttab), rownames(metatable))
# metatable = metatable[namecheck,]
# normcounttab = normcounttab[,namecheck]

data("sample_annotations")
data("sample_values")

comptable = metatable = sample_annotations
normcounttab = sample_values

outfilepath = "/Users/tosh/Desktop/Ruggles_Lab/projects/outlier-tool/output/"
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)


## EXTRACT COMPARISONS
groupings = comparison_groupings(comptable)

#####

## MAKE OUTLIER TABLE
#####
## Separate out the "i"th gene, take the bounds, and then create a column that says whether or not this gene is high, low, or none in a sample with regards to the other samples in the dataset. Repeat this for every gene to create a reference table
reftable_function_out = make_outlier_table(normcounttab,
                analyze_negative_outliers = TRUE, aggregate_features = TRUE,
                feature_delineator = "\\.")

outliertab = reftable_function_out$outliertab
upperboundtab = reftable_function_out$upperboundtab
lowerboundtab = reftable_function_out$lowerboundtab
sampmedtab = reftable_function_out$sampmedtab
aggnegoutlierstab = reftable_function_out$aggnegoutlierstab
aggnegfractiontab = reftable_function_out$aggnegfractiontab
aggposoutlierstab = reftable_function_out$aggposoutlierstab
aggposfractiontab = reftable_function_out$aggposfractiontab
#####

## Create Groupings and Test for Significance
## TABULATE OUTLIER DATA on a per gene basis for each of our subgroups
grouptablist1 = count_outliers(groupings, outliertab, analyze_negative_outliers = FALSE)
grouptablist2 = count_outliers(groupings, outliertab, analyze_negative_outliers = TRUE)
grouptablist3 = count_outliers(groupings, aggposoutlierstab, analyze_negative_outliers = FALSE)
grouptablist4 = count_outliers(groupings, aggnegoutlierstab, analyze_negative_outliers = TRUE)

## Create a blank list, and fill it with the pvalue for the fisher.test for each gene in the upward and downward significance directions
outlier_analysis_out1 = outlier_analysis(grouptablist1, write_out_tables = FALSE, outfilepath = outfilepath)
outlier_analysis_out2 = outlier_analysis(grouptablist2, write_out_tables = FALSE, outfilepath = outfilepath)
outlier_analysis_out3 = outlier_analysis(grouptablist3, write_out_tables = FALSE, outfilepath = outfilepath)
outlier_analysis_out4 = outlier_analysis(grouptablist4, write_out_tables = FALSE, outfilepath = outfilepath)

testnum_vec = c(1,2,3,4)
outlier_analysis_list = list(outlier_analysis_out1, outlier_analysis_out2, outlier_analysis_out3, outlier_analysis_out4)
for (analysiscount in seq_len(length(testnum_vec))){

    testnum = testnum_vec[analysiscount]
    outlier_analysis_param = outlier_analysis_list[[analysiscount]]

    for (tablecount in seq_len(length(outlier_analysis_param))) {
        write.table(outlier_analysis_param[[tablecount]],
                    paste0(outfilepath, names(outlier_analysis_param[tablecount]),
                           "_", testnum, ".csv"), row.names = FALSE, quote = FALSE, sep = ",")
    }
}


## For any significant genes, we then want to add a functionality that will plot a heatmap using the metadata and the genes significant via the outlier analysis - we want to run it for both directions
# analysis_num_param = 4
hm1 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out1, counttab = normcounttab,
                metatable = metatable, fdrcutoffvalue = 0.1, write_out_plot = FALSE, outfilepath = outfilepath, )
pdf(paste0(outfilepath, "test_hm1.pdf"))
hm1
junk<-dev.off()

hm2 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out2, counttab = normcounttab,
                metatable = metatable, fdrcutoffvalue = 0.1, write_out_plot = FALSE, outfilepath = outfilepath)
pdf(paste0(outfilepath, "test_hm2.pdf"))
hm2
junk<-dev.off()

hm3 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out3, counttab = aggposfractiontab,
                metatable = metatable, fdrcutoffvalue = 0.1, write_out_plot = FALSE, outfilepath = outfilepath)
pdf(paste0(outfilepath, "test_hm3.pdf"))
hm3
junk<-dev.off()

hm4 = outlier_heatmap(outlier_analysis_out = outlier_analysis_out4, counttab = aggnegfractiontab,
                      metatable = metatable, fdrcutoffvalue = 0.1, write_out_plot = FALSE, outfilepath = outfilepath)
pdf(paste0(outfilepath, "test_hm4.pdf"))
hm4
junk<-dev.off()

# for (i in 1:length(outlier_analysis_out)) {
#   analysis_num_param = i
#   metatablesort = metatable[order(metatable[,i]),]
#   outlier_heatmap(outlier_analysis_out, analysis_num = analysis_num_param, normcounttab, metatablesort, fdrcutoffvalue = 0.1)
# }


# temp1 = normcounttab[c("RECQL5-S727", "ARHGAP29-T948", "KMT2B-S844", "TMCC3-S253"),rownames(metatablesort)]
#
# maptab = as.matrix(t(apply(normcounttab[c("RECQL5-S727", "ARHGAP29-T948", "KMT2B-S844", "TMCC3-S253"),rownames(metatablesort)], 1, function(x) scale(x))))
# colnames(maptab) = colnames(normcounttab[c("RECQL5-S727", "ARHGAP29-T948", "KMT2B-S844", "TMCC3-S253"),rownames(metatablesort)])
# identical(colnames(maptab), colnames(normcounttab))



data("sample_values")
data("sample_annotations")

temp1 <- blacksheep(counttable = sample_values, metatable = sample_annotations,
    analyze_negative_outliers = TRUE, aggregate_features = TRUE,
    feature_delineator = "\\.", fdrcutoffvalue = 0.1,
    outfilepath = getwd())





