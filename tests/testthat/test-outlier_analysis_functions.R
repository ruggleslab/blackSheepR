context("blacksheepr Functions")
library(blacksheepr)

test_that("grouping function works", {
    data("sample_annotationdata")
    groupings = comparison_groupings(sample_annotationdata)

    expect_equal(sum(!is.na(sample_annotationdata)),
                 sum(sapply(groupings, length)))
})

test_that("making the outlier table works", {
    data("sample_phosphodata")
    reftable_function_out = make_outlier_table(sample_phosphodata[1:1000,])
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    sampmedtab = reftable_function_out$sampmedtab

    expect_equal(dim(outliertab), dim(sample_phosphodata[1:1000,]))
    expect_equal(sum(is.na(outliertab)),
                 sum(is.na(sample_phosphodata[1:1000,])))

    reftable_function_out = make_outlier_table(sample_phosphodata[1:1000,],
                                               analyze_negative_outliers = TRUE)
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    lowerboundtab = reftable_function_out$lowerboundtab
    sampmedtab = reftable_function_out$sampmedtab

    expect_true(!is.null(lowerboundtab))

})

test_that("counting outliers works", {
    data("sample_phosphodata")
    reftable_function_out = make_outlier_table(sample_phosphodata[1:1000,])
    outliertab = reftable_function_out$outliertab

    data("sample_annotationdata")
    groupings = comparison_groupings(sample_annotationdata)

    count_outliers_out = count_outliers(groupings, outliertab)
    grouptablist = count_outliers_out$grouptablist

    subcat = sort(as.vector(apply(sample_annotationdata, 2, function(x)
        na.omit(unique(x)))))

    expect_equal(length(grouptablist), length(subcat))

})


test_that("making the outlier table works", {
    data("sample_phosphodata")
    reftable_function_out = make_outlier_table(sample_phosphodata[1:1000,])
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    lowerboundtab = reftable_function_out$lowerboundtab
    sampmedtab = reftable_function_out$sampmedtab

    expect_equal(dim(outliertab), dim(sample_phosphodata[1:1000,]))
    expect_equal(sum(is.na(outliertab)),
                 sum(is.na(sample_phosphodata[1:1000,])))

})

test_that("outlier analysis works", {
    data("sample_phosphodata")
    reftable_function_out = make_outlier_table(sample_phosphodata[1:1000,])
    outliertab = reftable_function_out$outliertab

    data("sample_annotationdata")
    groupings = comparison_groupings(sample_annotationdata)

    grouptablist = count_outliers(groupings, outliertab)$grouptablist

    outlier_analysis_out = outlier_analysis(grouptablist)

    expect_equal(length(outlier_analysis_out), ncol(sample_annotationdata))

})


test_that("plotting the heatmap works", {
    data("sample_phosphodata")
    reftable_function_out = make_outlier_table(sample_phosphodata[1:1000,])
    outliertab = reftable_function_out$outliertab

    data("sample_annotationdata")
    groupings = comparison_groupings(sample_annotationdata)

    grouptablist = count_outliers(groupings, outliertab)$grouptablist

    outlier_analysis_out = outlier_analysis(grouptablist)

    hm1 = outlier_heatmap(outlier_analysis_out, analysis_num = NULL,
                    sample_phosphodata[1:1000,], sample_annotationdata,
                    fdrcutoffvalue = 0.1)

    ## Test to see if the number of analyses that have significant genes also
    ## have heatmaps
    expect_equal(length(hm1), sum(unlist(
        lapply(outlier_analysis_out, function(x)
            sum(rowSums(x[,c(4,5)] < 0.1) == 2) > 0))))

})
