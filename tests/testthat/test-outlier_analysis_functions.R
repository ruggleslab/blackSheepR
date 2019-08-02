context("BlackSheep Functions")
library(BlackSheep)

test_that("grouping function works", {
    data("sample_annotations")
    groupings = comparison_groupings(sample_annotations)

    expect_equal(sum(!is.na(sample_annotations)),
                 sum(sapply(groupings, length)))
})

test_that("making the outlier table works", {
    data("sample_values")
    reftable_function_out = make_outlier_table(sample_values)
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    sampmedtab = reftable_function_out$sampmedtab

    expect_equal(dim(outliertab), dim(sample_values))
    expect_equal(sum(is.na(outliertab)), sum(is.na(sample_values)))

    reftable_function_out = make_outlier_table(sample_values,
                                               analyze_negative_outliers = TRUE)
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    lowerboundtab = reftable_function_out$lowerboundtab
    sampmedtab = reftable_function_out$sampmedtab

    expect_true(!is.null(lowerboundtab))

    reftable_function_out = make_outlier_table(sample_values,
                                               analyze_negative_outliers = FALSE,
                                               aggregate_features = TRUE)
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    sampmedtab = reftable_function_out$sampmedtab
    aggposoutlierstab = reftable_function_out$aggposoutlierstab
    aggposfractiontab = reftable_function_out$aggposfractiontab

    expect_true(!is.null(aggposoutlierstab))
    expect_true(!is.null(aggposfractiontab))

    reftable_function_out = make_outlier_table(sample_values,
                                               analyze_negative_outliers = TRUE,
                                               aggregate_features = TRUE)
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    sampmedtab = reftable_function_out$sampmedtab
    lowerboundtab = reftable_function_out$lowerboundtab
    aggposoutlierstab = reftable_function_out$aggposoutlierstab
    aggposfractiontab = reftable_function_out$aggposfractiontab
    aggnegoutlierstab = reftable_function_out$aggnegoutlierstab
    aggnegfractiontab = reftable_function_out$aggnegfractiontab

    expect_true(!is.null(aggnegoutlierstab))
    expect_true(!is.null(aggnegfractiontab))

})

test_that("counting outliers works", {
    data("sample_values")
    reftable_function_out = make_outlier_table(sample_values)
    outliertab = reftable_function_out$outliertab

    data("sample_annotations")
    groupings = comparison_groupings(sample_annotations)

    grouptablist = count_outliers(groupings, outliertab)

    subcat = sort(as.vector(apply(sample_annotations, 2, function(x)
        na.omit(unique(x)))))

    expect_equal(length(grouptablist), length(subcat))

})


test_that("making the outlier table works", {
    data("sample_values")
    reftable_function_out = make_outlier_table(sample_values)
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    lowerboundtab = reftable_function_out$lowerboundtab
    sampmedtab = reftable_function_out$sampmedtab

    expect_equal(dim(outliertab), dim(sample_values))
    expect_equal(sum(is.na(outliertab)), sum(is.na(sample_values)))

})

test_that("outlier analysis works", {
    data("sample_values")
    reftable_function_out = make_outlier_table(sample_values)
    outliertab = reftable_function_out$outliertab

    data("sample_annotations")
    groupings = comparison_groupings(sample_annotations)

    grouptablist = count_outliers(groupings, outliertab)

    outlier_analysis_out = outlier_analysis(grouptablist)

    expect_equal(length(outlier_analysis_out), ncol(sample_annotations))

})


# test_that("plotting the heatmap works", {
#     data("sample_values")
#     reftable_function_out = make_outlier_table(sample_values)
#     outliertab = reftable_function_out$outliertab
#
#     data("sample_annotations")
#     groupings = comparison_groupings(sample_annotations)
#
#     grouptablist = count_outliers(groupings, outliertab)
#
#     outlier_analysis_out = outlier_analysis(grouptablist)
#
#     outlier_heatmap(outlier_analysis_out, analysis_num = NULL,
#                     sample_values, sample_annotations, fdrcutoffvalue = 0.1,
#                     outfilepath = getwd())
#
#
#
# })
