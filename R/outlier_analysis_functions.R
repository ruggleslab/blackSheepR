################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## blacksheep Script File


## MAKE COMPARISON COLUMNS
#' Utility function that will take in columns with several subcategories,
#' and output several columns each with binary classifications.
#' ex) col1: A,B,C >> colA: A,notA,notA; colB: notB,B,notB; colC: notC,notC,C
#'
#' @param intable table where each column has more than one subcategory, can
#'     be multiple columns
#' @return an expanded table with each of the column as a binary labeling of
#'     each subcatecory.
#' @keywords outliers
#' @import stats
#' @export
#' @examples
#' data("sample_annotations")
#' new_comparisons = make_comparison_columns(sample_annotations[,1,drop=FALSE])
make_comparison_columns <- function(intable){
    ## Create empty outtablelist
    outtablelist = list()
    ## Run through each inputted column
    for (i in seq_len(ncol(intable))) {
        ## Pull out the rownames and column names to save for later
        featurenames = rownames(intable)
        categoryname = colnames(intable)[i]

        ## Select the column we are working with
        subtable = intable[,i, drop=FALSE]
        subtable = apply(subtable, 2, as.character)

        ## Create duplicate columns - one for each subcategory
        new_comp_tab = matrix(rep(subtable, length(unique(subtable[,1]))),
                    nrow = nrow(intable), ncol = length(unique(subtable[,1])),
                    dimnames = list(featurenames, unique(subtable[,1])))

        ## For each subcategory, declare it as binary
        for ( ii in seq_len(ncol(new_comp_tab))) {
            new_comp_tab[new_comp_tab[,ii, drop=FALSE] !=
                            colnames(new_comp_tab[,ii, drop=FALSE]),ii] <-
                paste0("not_", colnames(new_comp_tab[,ii,drop=FALSE]))
        }

        ## Add back the column name, and then save out
        colnames(new_comp_tab) = paste0(categoryname, "_",
                                        colnames(new_comp_tab))
        outtablelist[[i]] = new_comp_tab

    }
    output = do.call(cbind, outtablelist)
    return(output)
}


## EXTRACT COMPARISONS
#' Create all of the groups based on the input metadata
#'
#' @param comptable table where each column will have comparisons drawn from it
#' @return a list with each of the groups as an entry in the list
#'     NOTE - this list will be ncol*2 long where ncol is the number comparisons
#' @keywords outliers
#' @import stats
#' @export
#' @examples
#' data("sample_annotations")
#' groupings = comparison_groupings(sample_annotations)
comparison_groupings <- function(comptable) {
    ## Initate blank list
    groupings = list()

    ## FAILSAFE - to insure that these are charactors not factors
    comptable = data.frame(comptable)
    samplenames = rownames(comptable)
    comptable = data.frame(lapply(comptable,as.character),
                                stringsAsFactors = FALSE)
    rownames(comptable) = samplenames
    for (i in seq_len(ncol(comptable))) {
        # Select column and create comparison based off entries in column
        subsamp = comptable[,i, drop=FALSE]
        subsamp[subsamp == ""] <- NA ## Failsafe to turn blank cells into NA
        subcats = unique(na.omit(subsamp[,1]))

        ## Fill the i'th entry in list with samples that fall into subcategory
        groupings[[i]] = rownames(na.omit(subsamp[subsamp == subcats[1],,
                                                drop=FALSE]))
        names(groupings)[i] = paste0(colnames(subsamp), "__", subcats[1])

        ## fill the i+num comps entry with samples in that subcategory
        groupings[[i+ncol(comptable)]] = rownames(na.omit(
            subsamp[subsamp == subcats[2],,drop=FALSE]))
        names(groupings)[i+ncol(comptable)] = paste0(
            colnames(subsamp), "__", subcats[2])
        ## Continue this through - this order enables a 2col matrix with these
        ## entries next to each other ex) 1,4;2,5;3,6 for 3 comps
}

    ## Define the output
    return(groupings)
}




## MAKE REF TABLE
#' Separate out the "i"th gene, take the bounds, and then create a column
#' that says whether or not this gene is high, low, or none in a sample with
#' regards to the other samples in the dataset. Repeat this for every gene to
#' create a reference table.
#' If Aggregating - this will output the aggregate count and fraction table
#' for your reference.
#'
#' @param intable table with all of the inputted information, samples along the
#'     x-axis, features along the y-axis
#' @param analyze_negative_outliers DEFAULT: FALSE; Toggle the analysis of
#'     outliers in the negative direction as well. Will lead to the output of
#'     the outlier table containing "-1" values, in addition to negative outputs
#'     for boundaries and aggregate tables (if applicable)
#' @return a list with varied sections depending on parameters:
#'     $outliertab - table converted to outlier form with 0s, 1s, and -1s,
#'     $upperboundtab - list of upper boundaries for outliers
#'     $lowerboundtab - list of lower boundaries of outliers
#'     $sampmedtab - list of median value per feature
#'
#' @keywords outliers
#' @export
#' @examples
#' data("sample_values")
#' reftable_function_out = make_outlier_table(sample_values,
#'     analyze_negative_outliers = TRUE)
#' outliertab = reftable_function_out$outliertab
#' upperboundtab = reftable_function_out$upperboundtab
#' lowerboundtab = reftable_function_out$lowerboundtab
#' sampmedtab = reftable_function_out$sampmedtab
make_outlier_table <- function(intable, analyze_negative_outliers = FALSE){

    # ## VECTORIZE THE PROCESS OF MAKING THE REFTABLE
    # medout = apply(intable, 1, function(x) median(x, na.rm = TRUE))
    # uppboundout = apply(intable, 1, function(x)
    #     median(x, na.rm = TRUE) + 1.5*IQR(x, na.rm = TRUE))
    # lowboundout = apply(intable, 1, function(x)
    #     median(x, na.rm = TRUE) - 1.5*IQR(x, na.rm = TRUE))
    #
    # set_pos_outliers <- function(x) {
    #     x[x > (median(x, na.rm = TRUE) + 1.5*IQR(x, na.rm = TRUE))] <- 1
    #     x[x < (median(x, na.rm = TRUE) + 1.5*IQR(x, na.rm = TRUE))] <- 0
    # }
    # outliertabout = do.call(rbind, lapply(intable, set_pos_outliers))

    ## Define blank out lists
    outlist = uppboundlist = lowboundlist = sampmedlist = list()
    for (intablegenecount in seq_len(nrow(intable))) {

        ## Save the sampname, then subset the table to just the one gene
        genename = rownames(intable)[intablegenecount]
        sampdat = as.data.frame(t(intable[intablegenecount, , drop=FALSE]))
        sampdat$categ = 0
        sampdat[is.na(sampdat[,1]),2] <- NA
        colnames(sampdat) = c(genename, genename)

        ## Set upper and lower bounds based off of the IQR and med
        ## Save upper and lower bound, and med to a list we then tabulate
        sampmedlist[[intablegenecount]] = sampmed =
            median(sampdat[,1],na.rm=TRUE)
        names(sampmedlist)[intablegenecount] = genename

        if (analyze_negative_outliers != TRUE) {
            uppboundlist[[intablegenecount]] = uppbound =
                sampmed + 1.5*IQR(sampdat[,1], na.rm = TRUE)
            names(sampmedlist)[intablegenecount] =
                names(uppboundlist)[intablegenecount] = genename
            sampdat[!is.na(sampdat[,1]) & sampdat[,1] > uppbound,2] <- 1
            lowboundlist = NULL
        } else {
            lowboundlist[[intablegenecount]] = lowbound =
                sampmed - 1.5*IQR(sampdat[,1], na.rm = TRUE)
            names(lowboundlist)[intablegenecount] = genename
            sampdat[!is.na(sampdat[,1]) & sampdat[,1] < lowbound,2] <- -1
            uppboundlist = NULL
        }


        ## Create a categ and assign high, low, none status given the boundaries
        # sampdat$categ = 0
        # sampdat[is.na(sampdat[,1]),2] <- NA
        # sampdat[!is.na(sampdat[,1]) & sampdat[,1] > uppbound,2] <- 1
        # colnames(sampdat) = c(genename, genename)

        ## Return negative information if param is TRUE
        # if (analyze_negative_outliers == TRUE) {
        #     # lowboundlist[[intablegenecount]] = lowbound =
        #     #     sampmed - 1.5*IQR(sampdat[,1], na.rm = TRUE)
        #     # names(lowboundlist)[intablegenecount] = genename
        #     # sampdat[!is.na(sampdat[,1]) & sampdat[,1] < lowbound,2] <- -1
        #     } else { lowboundlist = NULL}

        ## Save what we just made - which is n by 2 table with the gene in one
        ## col and the status of the sample in the other
        outlist[[intablegenecount]] = sampdat[,2, drop=FALSE]

    }
    ## Combine all of our lists
    # outliertab - the master table that tells whether each gene is significant
    # compared to the rest of the data set, upbound, downbound, sampmean - each
    # lists that have the upper boundary/lower boundary of signifigance/med for
    # the gene across all samples
    outliertab = as.data.frame(t(do.call(cbind, outlist)))
    # upperboundtab = do.call(rbind, uppboundlist)
    if (analyze_negative_outliers == TRUE) {
        negativeoutlierlist = list(lowerboundtab = do.call(rbind, lowboundlist))
        positiveoutlierlist = NULL
    } else {
        negativeoutlierlist = NULL
        positiveoutlierlist = list(upperboundtab = do.call(rbind, uppboundlist))
    }
    # lowerboundtab = do.call(rbind, lowboundlist)
    sampmedtab = do.call(rbind, sampmedlist)

    # if (analyze_negative_outliers == TRUE) {
    #     negativeoutlierlist = list(lowerboundtab = lowerboundtab)
    # } else {
    #     negativeoutlierlist = NULL
    # }
    output = c(list(
                outliertab = outliertab,
                sampmedtab = sampmedtab),
                positiveoutlierlist,
                negativeoutlierlist)

    ## Return the outputted values
    # outliertab where the intable has been turned into a reference table of
    # whether or not the gene is sig high, low, or none
    # upperboundtab - upperbound of signifigance for each of the genes
    # lowerboundtab - lowerbound of signifigance for each of the genes
    # sampmeantab - the average value for each of the genes
    return(output)
}


## TABULATE OUTLIER DATA on a per gene basis for each of our subgroups
#' Count up the outlier information for each of the groups you have made.
#' If aggregating then you will have to turn the parameter on, but you still
#'     input the outliertable. Aggregate will count the total number of
#'     outliers AND nonoutliers in its operation, so it needs the original
#'     outlier table made by the <make_outlier_table> function.
#' @usage count_outliers(groupings, outliertab,
#'     aggregate_features = FALSE, feature_delineator = "\\\\.")
#' @param groupings table generated by the comparison_groupings function
#' @param outliertab outlier table generated by make_outlier_table
#' @param aggregate_features DEFAULT: FALSE; Toggle the Aggregate feature, which
#'     will aggregate features in your table based on the given delineator.
#'     Aggregation will output counts for the TOTAL number of outliers and non-
#'     outliers across ALL sites you aggregate across.
#' @param feature_delineator DEFAULT: <"\\.">; What character delineates the
#'     separation between primary and secondary features. NOTE: to use proper
#'     R syntax with escape characters if necessary
#'     Ex) Protein1.Phosphosite1 uses "\\." to aggregate on Protein1
#' @return the tabulated information of outliers per group
#' @keywords outliers
#' @export
#' @examples
#'
#' data("sample_values")
#' reftable_function_out = make_outlier_table(sample_values)
#' outliertab = reftable_function_out$outliertab
#'
#' data("sample_annotations")
#' groupings = comparison_groupings(sample_annotations)
#'
#' count_outliers_out = count_outliers(groupings, outliertab,
#'     aggregate_features = FALSE)
#' grouptablist = count_outliers_out$grouptablist
#' fractiontab = count_outliers_out$fractiontab
count_outliers <- function(groupings, outliertab,
    aggregate_features = FALSE, feature_delineator = "\\.") {

    ## Define empty starting list
    grouptablist = list()
    ## Set factor depending on analysis - positive or negative
    #outliervalue = ifelse(analyze_negative_outliers == TRUE, -1, 1)

    ## AUTOMATICALLY DETECT OUTLIER VALUE
    outliervalue = unique(unlist(outliertab))[unique(unlist(outliertab)) != 0 &
                                    !is.na(unique(unlist(outliertab)))]
    ## AUTOMATICALLY DETECT OUTLIER VALUE

    if (aggregate_features == TRUE) {
        feature_labels = do.call(rbind, strsplit(sub(feature_delineator,
                                "xyz", rownames(outliertab)), split = "xyz"))
        colnames(feature_labels) = c("primary_feature", "secondary_feature")
        aggtabin = cbind(feature_labels, outliertab)

        aggoutliertab = aggregate(aggtabin[,c(3:ncol(aggtabin))],
                    by = list(primary_feature = aggtabin[,"primary_feature"]),
                    function(x) sum(x==outliervalue, na.rm = TRUE))
        aggoutliertab = data.frame(aggoutliertab[,-1],
                            row.names = aggoutliertab[,1], check.names = FALSE)

        aggnonoutliertab = aggregate(aggtabin[,c(3:ncol(aggtabin))],
                    by = list(primary_feature = aggtabin[,"primary_feature"]),
                    function(x) sum(x!=outliervalue, na.rm = TRUE))
        aggnonoutliertab = data.frame(aggnonoutliertab[,-1],
                        row.names = aggnonoutliertab[,1], check.names = FALSE)

        ## Create fraction table in similar manner
        fractiontab = aggregate(aggtabin[,c(3:ncol(aggtabin))],
            by = list(primary_feature = aggtabin[,"primary_feature"]),
                function(x) sum(x[x==outliervalue],na.rm = TRUE)/sum(!is.na(x)))
        fractiontab = data.frame(fractiontab[,-1],
                        row.names=fractiontab[,1], check.names = FALSE)
    } else {
        #if (analyze_negative_outliers == FALSE) {
        if (outliervalue == 1) {
            fractiontab = outliertab
        } else {
            fractiontab = -outliertab
        }
    }

    ## Perform counting in desired direction
    for (groupcount in seq_len(length(groupings))) {
        subgroup = outliertab[,colnames(outliertab) %in%
                                        groupings[[groupcount]], drop=FALSE]
        if (aggregate_features == FALSE) {
            grouptablist[[groupcount]] = list()
            grouptablist[[groupcount]][[1]] = t(apply(subgroup,MARGIN=1,
                        function(x) table(factor(x, levels=c(0,outliervalue)))))
            grouptablist[[groupcount]][[2]] = groupings[[groupcount]]
            names(grouptablist[[groupcount]]) = c("feature_counts", "samples")

        } else {
            aggoutliercount = cbind.data.frame(
                primary_feature = rownames(aggoutliertab),
                outlier_count = rowSums(
                                aggoutliertab[,groupings[[groupcount]]]))

            aggnonoutliercount = cbind.data.frame(
                primary_feature = rownames(aggnonoutliertab),
                outlier_count = rowSums(
                                aggnonoutliertab[,groupings[[groupcount]]]))

            aggcounttab = merge(aggnonoutliercount, aggoutliercount,
                                by = "primary_feature")
            rownames(aggcounttab) = aggcounttab[,1]
            colnames(aggcounttab) = c("primary_feature", 0, outliervalue)

            ## Save out the table to the first entry in the sublist
            grouptablist[[groupcount]] = list()
            grouptablist[[groupcount]][[1]] = aggcounttab[,c(2,3)]

            ## Also saving the samples as an additional entry in the sublist
            grouptablist[[groupcount]][[2]] = groupings[[groupcount]]
            names(grouptablist[[groupcount]]) = c("feature_counts", "samples")

        }
        names(grouptablist)[groupcount] = names(groupings)[groupcount]
    }
    ## Define our output with different lists depending on input params
    if (aggregate_features == TRUE){
        aggregateoutlist = list(aggoutliertab = aggoutliertab)
    } else {
        aggregateoutlist = NULL
    }
    output = c(list(grouptablist = grouptablist, fractiontab = fractiontab),
                    aggregateoutlist)
    return(output)
}


## DETERMINE P VALUE from the fisher.test for each gene in the upward and
## downward significance directions
#' With the grouptablist generated by count_outliers - run through and run a
#'     fisher exact test to get the p.value for the difference in outlier count
#'     for each feature in each of your comparisons
#' @usage outlier_analysis(grouptablist, fraction_table = NULL,
#'     fraction_samples_cutoff = 0.3,
#'     write_out_tables = FALSE, outfilepath = getwd())
#' @param grouptablist table generated by the count_outliers function. NOTE that
#'     the inputted grouptablist will be deciphered to determine its content.
#'     This means that user decides to input the outliertab or aggregate tab,
#'     and the output will analyze according to what positive and negative
#'     information is contained within the table
#' @param fraction_table DEFAULT: NULL; Input a fraction table to filter to
#'     only include features that have x% of samples in the ingroup that have
#'     an outlier.
#' @param fraction_samples_cutoff DEFAULT: 0.3; Input a fractional cut off for
#'     the of samples that need to have an outlier for feature to be
#'     considered. ex) 10 samples in ingroup - 3 need to have an outlier for
#'     feature to be considered significant
#' @param write_out_tables DEFAULT: FALSE; utility in function to write out
#'     each of the analyses to a separate table to whereever <outfilepath> is
#'     specfied.
#' @param outfilepath the full string path to where the file should output to,
#'     DEFAULT is current working directory
#' @return the analysis table with p.value, fdr, and raw data per comparison
#' @keywords outliers
#' @import stats utils
#' @export
#' @examples
#'
#' data("sample_values")
#' head(sample_values)
#' reftable_function_out = make_outlier_table(sample_values)
#' outliertab = reftable_function_out$outliertab
#'
#' data("sample_annotations")
#' groupings = comparison_groupings(sample_annotations)
#'
#' count_outliers_out = count_outliers(groupings, outliertab,
#'     aggregate_features = FALSE)
#' grouptablist = count_outliers_out$grouptablist
#' fractiontab = count_outliers_out$fractiontab
#'
#' outlier_analysis_out = outlier_analysis(grouptablist,
#'     fraction_table = fractiontab)
outlier_analysis <- function(grouptablist,
                        fraction_table = NULL, fraction_samples_cutoff = 0.3,
                        write_out_tables = FALSE, outfilepath = getwd()) {
    ## Define blank starting lists
    outtablelist = list()
    ## Create comparison matrix
    groupcombos = matrix(names(grouptablist),
                        nrow=(length(grouptablist)/2), ncol = 2)
    for (groupcombonum in seq_len(nrow(groupcombos))) {
        ## Pull out the two groups of interest for the run
        group1label = groupcombos[groupcombonum,1]
        group2label = groupcombos[groupcombonum,2]
        group1tab = as.data.frame(grouptablist[[group1label]][[1]])
        group2tab = as.data.frame(grouptablist[[group2label]][[1]])
        group1samps = grouptablist[[group1label]][[2]]
        group2samps = grouptablist[[group2label]][[2]]

        ## Run the function - but we want to run it 4 ways. Are there more pos
        ## outliers in A than B, are there more neg outliers in A than B, are
        ## there more pos outliers in B than A, are there more neg outliers in
        ## B than A

        ## Vectorized split_and_fish - extracting pval for each comparison
        comptab1 = merge(group2tab, group1tab, by = "row.names")
        comptab1a = data.frame(comptab1, row.names = comptab1[,1])[,2:5]
        comptab2 = merge(group1tab, group2tab, by = "row.names")
        comptab2a = data.frame(comptab2, row.names = comptab2[,1])[,2:5]

        split_and_fish <- function(combined_grouptab, sidedness){
            sigval = ifelse(sidedness=="up", 1, -1)
            conttab = matrix(unlist(combined_grouptab), nrow = 2, ncol = 2,
                        dimnames = list(c("0", sigval), c("group1","group2")))

            return(fisher.test(conttab, alternative = "two.sided")$p.value)
        }

        if ("1" %in% colnames(group1tab)) {
            upstatgroup1 = upstatgroup2 =
                apply(comptab1a, 1, split_and_fish, sidedness = "up")
        } else { upstatgroup1 = upstatgroup2 = NULL }

        if ("-1" %in% colnames(group1tab)) {
            downstatgroup1 = downstatgroup2 =
                apply(comptab1a, 1, split_and_fish, sidedness = "down")
        } else { downstatgroup1 = downstatgroup2 = NULL }

        fishout = cbind(upstatgroup1, upstatgroup2,
                        downstatgroup1, downstatgroup2)

        if (ncol(fishout) == 4) {colnames(fishout) = c(
            paste0("pval_more_pos_outliers_in_",
                    groupcombos[groupcombonum,c(1,2)]),
            paste0("pval_more_neg_outliers_in_",
                    groupcombos[groupcombonum,c(1,2)]))
        } else {
            if ("1" %in% colnames(group1tab)) {
                colnames(fishout)[c(1,2)] =
                    paste0("pval_more_pos_outliers_in_",
                            groupcombos[groupcombonum,c(1,2)])
            } else {colnames(fishout)[c(1,2)] =
                paste0("pval_more_neg_outliers_in_",
                        groupcombos[groupcombonum,c(1,2)]) }
        }

        #### RAW NUMBER FILTER
        ## Filter to only select features that already have a proportion of
        ## outliers greater in the ingroup
        group1prop_filter = (group1tab[,2]/(group1tab[,1] + group1tab[,2])) >
            (group2tab[,2]/(group2tab[,1] + group2tab[,2]))
        group2prop_filter = (group2tab[,2]/(group2tab[,1] + group2tab[,2])) >
            (group1tab[,2]/(group1tab[,1] + group1tab[,2]))
        group1prop_filter[is.na(group1prop_filter)] <- FALSE
        group2prop_filter[is.na(group2prop_filter)] <- FALSE

        ## Done separately so that if fraction table is incorporated,
        ## filters can be applied twice
        group1prop_filter_features = rownames(group1tab[group1prop_filter,])
        group2prop_filter_features = rownames(group2tab[group2prop_filter,])

        fishout = fishout[rownames(fishout) %in% union(
            group1prop_filter_features, group2prop_filter_features),,drop=FALSE]

        ## Fraction table filter
        if (!is.null(fraction_table)){
            group1fractab = fraction_table[,group1samps]
            group1fractab_select = group1fractab[rowSums(
                group1fractab!=0, na.rm = TRUE)/ncol(group1fractab) >
                    fraction_samples_cutoff,]
            group1fractab_select2 = group1fractab_select[
                rownames(group1fractab_select) %in% group1prop_filter_features,]

            group2fractab = fraction_table[,group2samps]
            group2fractab_select = group2fractab[rowSums(
                group2fractab!=0, na.rm = TRUE)/ncol(group2fractab) >
                    fraction_samples_cutoff,]
            group2fractab_select2 = group2fractab_select[
                rownames(group2fractab_select) %in% group2prop_filter_features,]

            fraction_selected_genes = union(rownames(group1fractab_select2),
                                            rownames(group2fractab_select2))
            fishout = fishout[rownames(fishout) %in% fraction_selected_genes,,
                            drop=FALSE]

        } else {fraction_selected_genes = rownames(fishout)}

        ## Add in FDR values for the pvalue metrics
        fish_to_fdr <- function(fishgroup){
            if (nrow(fishgroup) > 1) {fdrvals = apply(fishgroup, 2, function(x)
                    p.adjust(x, method = "BH"))
            } else {
                if(nrow(fishgroup) == 1){
                    fdrvals = fishgroup
                } else {
                    fdrvals = NULL
                }
            }
            return(fdrvals)
        }

        fdrgp1 = fishout[rownames(fishout) %in% group1prop_filter_features,
                        , drop = FALSE]
        fdrvals1 = fish_to_fdr(fdrgp1)
        fdrgp2 = fishout[rownames(fishout) %in% group2prop_filter_features,
                        , drop = FALSE]
        fdrvals2 = fish_to_fdr(fdrgp2)
        fdrvals = rbind(fdrvals1, fdrvals2)
        if (is.null(fdrvals)){next} # break to skip iteration if no sig genes
        colnames(fdrvals) = gsub(pattern = "pval", replacement = "fdr",
                                colnames(fdrvals))


    # ## Add in FDR values for the pvalue metrics
    # if (nrow(fishout) > 1) {
    #     ## APPLY FDR in UP and DOWN direction separately
    #     fdrgp1 = fishout[rownames(fishout) %in% group1prop_filter_features,]
    #     fdrvals1 = apply(fdrgp1, 2, function(x) p.adjust(x, method = "BH"))
    #     fdrgp2 = fishout[rownames(fishout) %in% group2prop_filter_features,]
    #     fdrvals2 = apply(fdrgp2, 2, function(x) p.adjust(x, method = "BH"))
    #     fdrvals = rbind(fdrvals1, fdrvals2)
    #
    #     colnames(fdrvals) = gsub(pattern = "pval", replacement = "fdr",
    #                              colnames(fdrvals))
    # } else {
    #     if (nrow(fishout) == 1){
    #         fdrvals = data.frame(as.list(p.adjust(fishout, method = "BH")),
    #                             row.names = rownames(fishout))
    #         colnames(fdrvals) = gsub(pattern = "pval", replacement = "fdr",
    #                                  colnames(fishout))
    #     } else {
    #         next
    #     }
    # }

        ## SAVE OUT DATA
        if ("1" %in% colnames(group1tab)) {
            colnames(group1tab) = paste(group1label, "_",
                                    c("nonposoutlier","posoutlier"), sep="")
            colnames(group2tab) = paste(group2label, "_",
                                    c("nonposoutlier","posoutlier"), sep="")
        } else{
            colnames(group1tab) = paste(group1label, "_",
                                    c("nonnegoutlier","negoutlier"), sep="")
            colnames(group2tab) = paste(group2label, "_",
                                    c("nonnegoutlier","negoutlier"), sep="")
        }

        ## Write out the final table
        outtablefin = Reduce(function(dtf1, dtf2) merge(dtf1, dtf2,by = "gene"),
                        lapply(list(fishout, fdrvals,
                        group1tab[fraction_selected_genes,],
                        group2tab[fraction_selected_genes,]),
                        function(x) data.frame(x, gene = row.names(x))))

        ## Format the final table to only have the pvalues for the analysis that
        ## leads to more in ingroup/outgroup
        outtablefin[outtablefin[,7]/(outtablefin[,6] + outtablefin[,7]) >
                outtablefin[,9]/(outtablefin[,8] + outtablefin[,9]),c(3,5)] = ""
        outtablefin[outtablefin[,7]/(outtablefin[,6] + outtablefin[,7]) <
                outtablefin[,9]/(outtablefin[,8] + outtablefin[,9]),c(2,4)] = ""

        ## Write out the tables
        if (write_out_tables == TRUE) {
            outlieranalysisoutfile = paste(outfilepath, "outlieranalysis_for_",
                group1label, "_vs_",
                group2label, ".csv", sep="")
            write.table(outtablefin, outlieranalysisoutfile, quote = FALSE,
                        sep = ",", row.names = FALSE)
        }

        ## Save out results in a list
        outtablelist[[groupcombonum]] = outtablefin
        names(outtablelist)[groupcombonum] = paste("outlieranalysis_for_",
            group1label, "_vs_",
            group2label, sep="")
    }
    return(Filter(Negate(is.null), outtablelist))
}




## PLOT HEATMAP with metadata, original countdata, and the outlieranalysis to
## pull out significant genes
#' With the grouptablist generated by count_outliers - run through and run a
#'     fisher exact test to get the p.value for the difference in outlier count
#'     for each feature in each of your comparisons
#' @usage outlier_heatmap(outlier_analysis_out, analysis_num = NULL,
#'     counttab, metatable, fdrcutoffvalue = 0.1,
#'     write_out_plot = FALSE, outfilepath = getwd())
#' @param outlier_analysis_out the full outlier_analysis data objet
#' @param analysis_num DEFAULT: NULL; if you only want to plot the heatmap for
#'     a particular analysis, enter number of that analysis
#' @param counttab the raw data before outlier analysis
#' @param metatable the complete metatable that was used to generate the
#'     comparisons, will be used for annotation of the heatmap
#' @param fdrcutoffvalue DEFAULT: 0.1; The FDR value for significance
#' @param write_out_plot DEFAULT: FALSE; write out the plot to <outfilepath>
#' @param outfilepath the full string path to where the file should output to,
#'     DEFAULT is current working directory
#' @return outputs a pdf with the heatmap in the current working directory
#' @keywords outliers
#' @import ComplexHeatmap RColorBrewer circlize
#' @export
#' @examples
#'
#' data("sample_values")
#' reftable_function_out = make_outlier_table(sample_values)
#' outliertab = reftable_function_out$outliertab
#'
#' data("sample_annotations")
#' groupings = comparison_groupings(sample_annotations)
#'
#' count_outliers_out = count_outliers(groupings, outliertab,
#'     aggregate_features = FALSE)
#' grouptablist = count_outliers_out$grouptablist
#' fractiontab = count_outliers_out$fractiontab
#'
#' outlier_analysis_out = outlier_analysis(grouptablist,
#'     fraction_table = fractiontab)
#'
#' metatable = sample_annotations
#' counttab = sample_values
#'
#' hm1 = outlier_heatmap(outlier_analysis_out, analysis_num = NULL,
#'     fractiontab, metatable, fdrcutoffvalue = 0.1)
outlier_heatmap <- function(outlier_analysis_out, analysis_num = NULL, counttab,
                            metatable, fdrcutoffvalue = 0.1,
                            write_out_plot = FALSE, outfilepath = getwd()) {
    ## Pull out the comparison columns from the metadata
    #compcols = colnames(metatable)[grep("comp_", colnames(metatable))]

    ## Define the start and stops of the foorloop - if they put in a single
    ## analysis to do - it will only output that, otherwise, will loop over all
    startcount = ifelse(is.null(analysis_num), 1, analysis_num)
    endcount = ifelse(is.null(analysis_num), length(outlier_analysis_out),
                    analysis_num)
    heatmaplist = NULL
    for (analysiscount in startcount:endcount) {
        intable = outlier_analysis_out[[analysiscount]]

        ## Select the column from the outlier_analysis_out that contain "fdr"
        ## and then grab rownames for columns that have sig value
        fdrcols = intable[ ,grepl( "fdr", colnames(intable))]
        GOI = intable[rowSums(fdrcols < fdrcutoffvalue) == 2,1]

        if (length(GOI) > 0) {
            ## Take the metatable, order it by 1s and then 2s on whatever
            ## comparison we are doing, take comparison columns for plotting
            ## Added in as.character as a failsafe
            subsetcounttab = counttab[as.character(GOI),
                                        rownames(metatable), drop=FALSE]

            annotation1 = annotationlist_builder(metatable)
            outfile1 = paste(outfilepath, "outlieranalysis_heatmap_for_",
                names(outlier_analysis_out)[analysiscount], ".pdf", sep="")
            hm1 = create_heatmap(counttab = subsetcounttab,
                colmetatable = metatable, colannotationlist = annotation1,
                colclusterparam = FALSE, rowclusterparam = FALSE,
                write_out_plot = write_out_plot,
                nameparam = names(outlier_analysis_out)[analysiscount],
                pdfoutfile = outfile1)
            #return(hm1)
            heatmaplist[analysiscount] = list(hm1)
            names(heatmaplist)[analysiscount] = paste0("print_",
                                names(outlier_analysis_out)[analysiscount])
        } else {
            #print(paste("No Significant Outliers for ",
            #names(outlier_analysis_out)[analysiscount],
            #" at an FDR cut off value of ", fdrcutoffvalue, sep = ""))
            }
    }
    return(Filter(Negate(is.null), heatmaplist))
}


## COMPLETE BLACKSHEEP FUNCTION
#' Run the entire blacksheep Function from Start to finish
#' @usage deva(se, analyze_negative_outliers = FALSE,
#'     aggregate_features = FALSE, feature_delineator = "\\\\.",
#'     fraction_samples_cutoff = 0.3, fdrcutoffvalue = 0.01,
#'     write_out = FALSE, outfilepath = getwd())
#' @param se The SummarizedExperiment object containing the countdata and the
#'     associated annotation data with comparisons in the colData object.
#' @param analyze_negative_outliers DEFAULT: FALSE; Toggle the analysis of
#'     outliers in the negative direction as well. Will lead to the output of
#'     the outlier table containing "-1" values, in addition to negative outputs
#'     for boundaries and aggregate tables (if applicable)
#' @param aggregate_features DEFAULT: FALSE; Toggle the Aggregate feature, which
#'     will aggregate features in your table based on the given delineator.
#'     Aggregation will output an aggregate table that counts the number of
#'     outliers per feature, and also a fraction table that show the number of
#'     outliers / number of candidates (which excludes missing values)
#' @param feature_delineator DEFAULT: <"\\."> What character delineates the
#'     separation between primary and secondary features. NOTE: to use proper
#'     R syntax with escape characters if necessary
#'     Ex) Protein1.Phosphosite1 uses "\\." to aggregate on Protein1
#' @param fraction_samples_cutoff DEFAULT: 0.3; Input a fractional cut off for
#'     the of samples that need to have an outlier for feature to be
#'     considered. ex) 10 samples in ingroup - 3 need to have an outlier for
#'     feature to be considered significant
#' @param fdrcutoffvalue DEFAULT: 0.1; The FDR value for significance
#' @param write_out DEFAULT: FALSE; write out tables and plots
#' @param outfilepath the full string path to where the file should output to,
#'     DEFAULT is current working directory
#' @return outputs a pdf with the heatmap in the current working directory
#' @keywords outliers
#' @import ComplexHeatmap RColorBrewer circlize
#' @rawNamespace import(SummarizedExperiment, except = c(start, end))
#' @export
#' @examples
#'
#' library(SummarizedExperiment)
#' data("sample_values")
#' data("sample_annotations")
#'
#' se = SummarizedExperiment(assays = list(counts = as.matrix(sample_values)),
#'     colData = DataFrame(sample_annotations))
#'
#' deva(se = se,
#'     analyze_negative_outliers = TRUE, aggregate_features = TRUE,
#'     feature_delineator = "\\.", fdrcutoffvalue = 0.1, write_out = FALSE,
#'     outfilepath = getwd())
deva <- function(se, analyze_negative_outliers = FALSE,
                        aggregate_features = FALSE, feature_delineator = "\\.",
                        fraction_samples_cutoff = 0.3, fdrcutoffvalue = 0.01,
                        write_out = FALSE, outfilepath = getwd()) {

    ## Extracting information from the SummarizedExperiment
    counttable = assays(se)[[1]]
    metatable = colData(se)

    ## Use the groupings function to create comparison groups
    groupings = comparison_groupings(metatable)

    ## Create the outlier table and other outlier objects
    reftable_function_out = make_outlier_table(counttable,
                analyze_negative_outliers = analyze_negative_outliers)
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    lowerboundtab = reftable_function_out$lowerboundtab
    sampmedtab = reftable_function_out$sampmedtab

    ## If there is no aggregation of features - then just run through normal
    ## analysis, outputting whats parameterized, withoutput as hm and tables
    ## Running through positive and negative side if necessary
    if (aggregate_features == FALSE) {
        ## Positive/nonaggregated workflow
        pos_count_outliers_out = count_outliers(groupings, outliertab)
        posgrouptablist = pos_count_outliers_out$grouptablist
        posfractiontab = pos_count_outliers_out$fractiontab

        pos_outlier_analysis_out = outlier_analysis(
            grouptablist = posgrouptablist,
            fraction_table = posfractiontab, fraction_samples_cutoff = 0.3,
            write_out_tables = write_out, outfilepath = outfilepath)

        hm1 = outlier_heatmap(outlier_analysis_out = pos_outlier_analysis_out,
                        analysis_num = NULL, counttab = posfractiontab,
                        metatable = metatable, fdrcutoffvalue = fdrcutoffvalue,
                        write_out_plot = write_out, outfilepath = outfilepath)

        ## Expanded out Code to for loop to sort the metatable for each comp
        hm1list = list()
        for (i in seq_len(length(pos_outlier_analysis_out))){
            plottable = metatable[do.call(order, c(decreasing = TRUE,
                                data.frame(metatable[,c(i,
                                    setdiff(seq_len(ncol(metatable)),i))]))),]
            hm1list[[i]] = outlier_heatmap(
                outlier_analysis_out = pos_outlier_analysis_out,
                analysis_num = i, counttab = posfractiontab,
                metatable = plottable, fdrcutoffvalue = fdrcutoffvalue,
                write_out_plot = write_out, outfilepath = outfilepath)
        }
        hm1 = unlist(hm1list)

        ## Negative/nonaggregated workflow
        if (analyze_negative_outliers == TRUE) {
            neg_count_outliers_out = count_outliers(groupings, outliertab,
                aggregate_features = FALSE)
            neggrouptablist = neg_count_outliers_out$grouptablist
            negfractiontab = neg_count_outliers_out$fractiontab

            neg_outlier_analysis_out = outlier_analysis(
                grouptablist = neggrouptablist,
                fraction_table = negfractiontab, fraction_samples_cutoff = 0.3,
                write_out_tables = write_out, outfilepath = outfilepath)

            ## Expanded out Code to for loop to sort the metatable for each comp
            hm2list = list()
            for (i in seq_len(length(pos_outlier_analysis_out))){
                plottable = metatable[do.call(order, c(decreasing = TRUE,
                                    data.frame(metatable[,c(i,
                                    setdiff(seq_len(ncol(metatable)),i))]))),]
                hm2list[[i]] = outlier_heatmap(
                    outlier_analysis_out = neg_outlier_analysis_out,
                    analysis_num = i, counttab = negfractiontab,
                    metatable = plottable, fdrcutoffvalue = fdrcutoffvalue,
                    write_out_plot = write_out, outfilepath = outfilepath)
            }
            hm2 = unlist(hm2list)
        }

        ## Return the output - parameterized to output pos/neg as appropriate
        if (analyze_negative_outliers != TRUE) {
            return(list(pos_outlier_analysis = pos_outlier_analysis_out,
                        significant_pos_heatmaps = hm1))
        } else {
            return(list(pos_outlier_analysis = pos_outlier_analysis_out,
                        significant_pos_heatmaps = hm1,
                        neg_outlier_analysis = neg_outlier_analysis_out,
                        significant_neg_heatmaps = hm2))
        }

    ## With feature aggregation - we do positive and negative separtely
    } else {
        pos_count_outliers_out = count_outliers(groupings, outliertab,
            aggregate_features = TRUE, feature_delineator = feature_delineator)
        posgrouptablist = pos_count_outliers_out$grouptablist
        posfractiontab = pos_count_outliers_out$fractiontab

        pos_outlier_analysis_out = outlier_analysis(
            grouptablist = posgrouptablist,
            fraction_table = posfractiontab,
            fraction_samples_cutoff = fraction_samples_cutoff,
            write_out_tables = write_out, outfilepath = outfilepath)

        ## Expanded out Code to for loop to sort the metatable for each comp
        hm1list = list()
        for (i in seq_len(length(pos_outlier_analysis_out))){
            plottable = metatable[do.call(order, c(decreasing = TRUE,
                                data.frame(metatable[,c(i,
                                setdiff(seq_len(ncol(metatable)),i))]))),]
            hm1list[[i]] = outlier_heatmap(
                outlier_analysis_out = pos_outlier_analysis_out,
                analysis_num = i, counttab = posfractiontab,
                metatable = plottable, fdrcutoffvalue = fdrcutoffvalue,
                write_out_plot = write_out, outfilepath = outfilepath)
        }
        hm1 = unlist(hm1list)

        if (analyze_negative_outliers == TRUE) {
            neg_count_outliers_out = count_outliers(groupings, outliertab,
                aggregate_features = TRUE,
                feature_delineator = feature_delineator)
            neggrouptablist = neg_count_outliers_out$grouptablist
            negfractiontab = neg_count_outliers_out$fractiontab

            neg_outlier_analysis_out = outlier_analysis(
                grouptablist = neggrouptablist,
                fraction_table = negfractiontab,
                fraction_samples_cutoff = fraction_samples_cutoff,
                write_out_tables = write_out, outfilepath = outfilepath)

            ## Expanded out Code to for loop to sort the metatable for each comp
            hm2list = list()
            for (i in seq_len(length(pos_outlier_analysis_out))){
                plottable = metatable[do.call(order, c(decreasing = TRUE,
                                    data.frame(metatable[,c(i,
                                    setdiff(seq_len(ncol(metatable)),i))]))),]
                hm2list[[i]] = outlier_heatmap(
                    outlier_analysis_out = neg_outlier_analysis_out,
                    analysis_num = i, counttab = negfractiontab,
                    metatable = plottable, fdrcutoffvalue = fdrcutoffvalue,
                    write_out_plot = write_out, outfilepath = outfilepath)
            }
            hm2 = unlist(hm2list)
        }

        if (analyze_negative_outliers != TRUE) {
            return(list(pos_outlier_analysis = pos_outlier_analysis_out,
                        significant_pos_heatmaps = hm1))
        } else {
            return(list(pos_outlier_analysis = pos_outlier_analysis_out,
                        significant_pos_heatmaps = hm1,
                        neg_outlier_analysis = neg_outlier_analysis_out,
                        significant_neg_heatmaps = hm2))
        }
    }

}

