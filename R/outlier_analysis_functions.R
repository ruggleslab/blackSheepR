################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Blacksheep Script File


## EXTRACT COMPARISONS
#####
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
    groupings = list()
    for (i in seq_len(ncol(comptable))) {
        # Select column and create coparison based off entries in column
        subsamp = comptable[,i, drop=FALSE]
        subsamp[subsamp == ""] <- NA ## Failsafe to turn blank cells into NA
        subcats = unique(na.omit(subsamp[,1]))

        groupings[[i]] = rownames(na.omit(subsamp[
            subsamp == subcats[1],,drop=FALSE]))
        names(groupings)[i] = paste0(colnames(subsamp), "__", subcats[1])

        groupings[[i+ncol(comptable)]] = rownames(na.omit(
            subsamp[subsamp == subcats[2],,drop=FALSE]))
        names(groupings)[i+ncol(comptable)] = paste0(
            colnames(subsamp), "__", subcats[2])
}

    ## Define the output
    return(groupings)
}
#groupings = comparison_groupings(comptable)
#####


## MAKE REF TABLE
#' Separate out the "i"th gene, take the bounds, and then create a column
#' that says whether or not this gene is high, low, or none in a sample with
#' regards to the other samples in the dataset. Repeat this for every gene to
#' create a reference table
#'
#' @param intable table with all of the inputted information, samples along the
#'     x-axis, features along the y-axis
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
#' @return a list with varied sections depending on parameters:
#'     $outliertab - table converted to outlier form with 0s, 1s, and -1s,
#'     $upperboundtab - list of upper boundaries for outliers
#'     $lowerboundtab - list of lower boundaries of outliers
#'     $sampmedtab - list of median value per feature
#'     $aggposoutlierstab - aggregate table for positive outliers
#'     $aggnegoutlierstab - aggregate table for negative outliers
#'     $aggposfractiontab - fraction table for number of pos outliers/candidates
#'     $aggnegfractiontab - fraction table for number of neg outliers/candidates
#'
#' @keywords outliers
#' @export
#' @examples
#' data("sample_values")
#' reftable_function_out = make_outlier_table(sample_values,
#'     analyze_negative_outliers = TRUE, aggregate_features = TRUE,
#'     feature_delineator = "\\.")
#' outliertab = reftable_function_out$outliertab
#' upperboundtab = reftable_function_out$upperboundtab
#' lowerboundtab = reftable_function_out$lowerboundtab
#' sampmedtab = reftable_function_out$sampmedtab
#' aggposoutlierstab = reftable_function_out$aggposoutlierstab
#' aggnegoutlierstab = reftable_function_out$aggnegoutlierstab
#' aggposfractiontab = reftable_function_out$aggposfractiontab
#' aggnegfractiontab = reftable_function_out$aggnegfractiontab
make_outlier_table <- function(intable, analyze_negative_outliers = FALSE,
                            aggregate_features = FALSE,
                            feature_delineator = "\\.") {

    outlist = uppboundlist = lowboundlist = sampmedlist = list()
    for (intablegenecount in seq_len(nrow(intable))) {

        ## Save the sampname, then subset the table to just the one sample
        genename = rownames(intable)[intablegenecount]
        sampdat = as.data.frame(t(intable[intablegenecount, , drop=FALSE]))

        ## Set upper and lower bounds based off of the IQR and med
        ## Save upper and lower bound, and med to a list we then tabulate
        sampmedlist[[intablegenecount]] = sampmed = median(sampdat[,1],
                                                            na.rm=TRUE)
        uppboundlist[[intablegenecount]] = uppbound = sampmed +
            1.5*IQR(sampdat[,1], na.rm = TRUE)
        names(sampmedlist)[intablegenecount] =
            names(uppboundlist)[intablegenecount] = genename

        ## Create a categ and assign high, low, none status given the boundaries
        sampdat$categ = 0
        sampdat[is.na(sampdat[,1]),2] <- NA
        sampdat[!is.na(sampdat[,1]) & sampdat[,1] > uppbound,2] <- 1
        colnames(sampdat) = c(genename, genename)

        if (analyze_negative_outliers == TRUE) {
            lowboundlist[[intablegenecount]] = lowbound = sampmed -
                1.5*IQR(sampdat[,1], na.rm = TRUE)
            names(lowboundlist)[intablegenecount] = genename
            sampdat[!is.na(sampdat[,1]) & sampdat[,1] < lowbound,2] <- -1
            } else { lowboundlist = NULL}

        ## Save what we just made - which is n by 2 table with the gene in one
        ## col and the status of the sample in the other
        outlist[[intablegenecount]] = sampdat[,2, drop=FALSE]

    }
    ## Combine all of our lists
    # reftable - the master table that tells whether each gene is significant
    # compared to the rest of the data set, upbound, downbound, sampmean - each
    # lists that have the upper boundary/lower boundary of signifigance/med for
    # the gene across all samples
    outliertab = as.data.frame(t(do.call(cbind, outlist)))
    upperboundtab = do.call(rbind, uppboundlist)
    if (analyze_negative_outliers == TRUE) {
        lowerboundtab = do.call(rbind, lowboundlist)}
    sampmedtab = do.call(rbind, sampmedlist)

    ## ADD IN AGGREGATE FUNCTION
    if (aggregate_features == TRUE) {
        # Set the delineator for feature split
        feature_delineator = "\\."
        # Reassign the FIRST instance of the delineator with a placeholder
        # Then split on that placeholder, and rename the columns to attach
        feature_labels = do.call(rbind, strsplit(sub(feature_delineator,
                        "xyz123", rownames(outliertab)), split = "xyz123"))
        colnames(feature_labels) = c("primary_feature", "secondary_feature")
        aggtabin = cbind(feature_labels, outliertab)

        tempfunc <- function(x) {
            print(x)
            sum(x[x==1], na.rm=TRUE)
        }


        make_aggregate_tables<- function(aggtabin, outliervalue) {
            aggoutliertab = aggregate(aggtabin[,c(3:ncol(aggtabin))],
                    by = list(primary_feature = aggtabin[,"primary_feature"]),
                    function(x) sum(x[x==outliervalue], na.rm = TRUE))
            aggoutliertab = data.frame(aggoutliertab[,-1],
                                        row.names=aggoutliertab[,1])

            aggfractiontab = aggregate(aggtabin[,c(3:ncol(aggtabin))],
                by = list(primary_feature = aggtabin[,"primary_feature"]),
                function(x)
                    sum(x[x==outliervalue], na.rm = TRUE)/sum(!is.na(x)))
            aggfractiontab = data.frame(aggfractiontab[,-1],
                                        row.names=aggfractiontab[,1])
            return(list(aggoutliertab = aggoutliertab,
                        aggfractiontab = aggfractiontab))
        }
        posaggresults = make_aggregate_tables(aggtabin, 1)
        if (analyze_negative_outliers == TRUE) {
            negaggresults = make_aggregate_tables(aggtabin, -1)
        }
    }

    if (aggregate_features == TRUE){
        aggregatelist = list(aggposoutlierstab = posaggresults$aggoutliertab,
                            aggposfractiontab = posaggresults$aggfractiontab)
    } else {
        aggregatelist = NULL
    }
    if (analyze_negative_outliers == TRUE) {
        negativeoutlierlist = list(lowerboundtab = lowerboundtab)
        if (aggregate_features == TRUE) {negativeoutlierlist = list(
            lowerboundtab = lowerboundtab,
            aggnegoutlierstab = negaggresults$aggoutliertab,
            aggnegfractiontab = negaggresults$aggfractiontab)
        }
    } else {
        negativeoutlierlist = NULL
    }
    output = c(list(
                outliertab = outliertab,
                sampmedtab = sampmedtab,
                upperboundtab = upperboundtab),
                negativeoutlierlist,
                aggregatelist
                )

    ## Return the outputted values
    # outliertab where the intable has been turned into a reference table of
    # whether or not the gene is sig high, low, or none
    # upboundtab - upperbound of signifigance for each of the genes
    # downboundtab - lowerbound of signifigance for each of the genes
    # sampmeantab - the average value for each of the genes
    return(output)
}


## TABULATE OUTLIER DATA on a per gene basis for each of our subgroups
#' Count up the outlier information for each of the groups you have made
#' @param groupings table generated by the comparison_groupings function
#' @param outliertab outlier table generated by make_outlier_table
#' @param analyze_negative_outliers DEFAULT: FALSE; analyze negative outliers in
#'     addition to positive outliers, NOTE - this MUST BE SET TO TRUE if using
#'     an aggregate negative table to get a result
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
#' grouptablist = count_outliers(groupings, outliertab,
#'     analyze_negative_outliers = FALSE)
count_outliers <- function(groupings, outliertab,
                            analyze_negative_outliers = FALSE) {
    grouptablist = list()
    if(analyze_negative_outliers == TRUE) {levelsparam = c(0,-1)} else
        {levelsparam = c(0,1)}
    for (groupcount in seq_len(length(groupings))) {
        subgroup = outliertab[,colnames(outliertab) %in%
                                        groupings[[groupcount]], drop=FALSE]
        grouptablist[[groupcount]] = t(apply(subgroup,MARGIN=1, function(x)
            table(factor(x, levels=levelsparam))))
        names(grouptablist)[groupcount] = names(groupings)[groupcount]
    }
    return(grouptablist)
}


## DETERMINE P VALUE from the fisher.test for each gene in the upward and
## downward significance directions
#' With the grouptablist generated by count_outliers - run through and run a
#'     fisher exact test to get the p.value for the difference in outlier count
#'     for each feature in each of your comparisons
#' @param grouptablist table generated by the count_outliers function. NOTE that
#'     the inputted grouptablist will be deciphered to determine its content.
#'     This means that user decides to input the outliertab or aggregate tab,
#'     and the output will analyze according to what positive and negative
#'     information is contained within the table
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
#' grouptablist = count_outliers(groupings, outliertab)
#'
#' outlier_analysis_out = outlier_analysis(grouptablist)
outlier_analysis <- function(grouptablist, write_out_tables = FALSE,
                            outfilepath = getwd()) {
    outtablelist = statoutlist = list()
    groupcombos = matrix(names(grouptablist),
                        nrow=(length(grouptablist)/2), ncol = 2)
    for (groupcombonum in seq_len(nrow(groupcombos))) {
        ## Pull out the two groups of interest for the run
        group1tab = as.data.frame(
            grouptablist[groupcombos[groupcombonum,1]][[1]])
        group2tab = as.data.frame(
            grouptablist[groupcombos[groupcombonum,2]][[1]])

        ## Run the function - but we want to run it 4 ways. Are there more pos
        ## outliers in A than B, are there more neg outliers in A than B, are
        ## there more pos outliers in B than A, are there more neg outliers in
        ## B than A
        for (fishtestgene in seq_len(nrow(group1tab))) {
            # function that will do the work of dividing the subgroup tables
            # into up and down directions, and perform the fisher test group1tab
            # and group2tab are binary tables for the subgroups, count is the
            # gene count (to work with the loop) and sidedness if "up"/"down"
            split_and_fish <- function(group1tab, group2tab, count, sidedness){
                if (sidedness=="up"){sigval = "1"}
                if (sidedness=="down"){sigval = "-1"}
                conttab = cbind(t(group1tab[count,c(sigval,"0"), drop=FALSE]),
                                t(group2tab[count,c(sigval,"0"), drop=FALSE]))
                colnames(conttab) = c("group1","group2")
                return(fisher.test(conttab, alternative = "greater")$p.value)
            }

            if ("1" %in% colnames(group1tab)) {
                upstatgroup1 = split_and_fish(group1tab, group2tab,
                                            fishtestgene, "up")
                upstatgroup2 = split_and_fish(group2tab, group1tab,
                                            fishtestgene, "up")
            } else { upstatgroup1 = upstatgroup2 = NULL }

            if ("-1" %in% colnames(group1tab)) {
                downstatgroup1 = split_and_fish(
                    group1tab, group2tab, fishtestgene, "down")
                downstatgroup2 = split_and_fish(
                    group2tab, group1tab, fishtestgene, "down")
            } else { downstatgroup1 = downstatgroup2 = NULL }
            statoutlist[[fishtestgene]] = c(upstatgroup1, upstatgroup2,
                                            downstatgroup1, downstatgroup2)
            names(statoutlist)[fishtestgene] = rownames(group1tab)[fishtestgene]
        }

        ## Write out the values from the fisher test
        fishout = do.call(rbind, statoutlist)

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

        ## Add in FDR values for the pvalue metrics
        fdrvals = apply(fishout, 2, function(x) p.adjust(x, method = "BH"))
        colnames(fdrvals) = gsub(pattern = "pval", replacement = "fdr",
                                colnames(fdrvals))


        ## SAVE OUT DATA
        if (ncol(group1tab) == 3){
            colnames(group1tab) = paste(groupcombos[groupcombonum,1], "_",
                            c("negoutlier","nonoutlier","posoutlier"), sep="")
            colnames(group2tab) = paste(groupcombos[groupcombonum,2], "_",
                            c("negoutlier","nonoutlier","posoutlier"), sep="")
        } else {
            if ("1" %in% colnames(group1tab)) {
                colnames(group1tab) = paste(groupcombos[groupcombonum,1], "_",
                                        c("nonposoutlier","posoutlier"), sep="")
                colnames(group2tab) = paste(groupcombos[groupcombonum,2], "_",
                                        c("nonposoutlier","posoutlier"), sep="")
            } else{
                colnames(group1tab) = paste(groupcombos[groupcombonum,1], "_",
                                        c("nonnegoutlier","negoutlier"), sep="")
                colnames(group2tab) = paste(groupcombos[groupcombonum,2], "_",
                                        c("nonnegoutlier","negoutlier"), sep="")
            }
        }
        outtablefin = cbind(gene = rownames(fishout),fishout, fdrvals,
                            group1tab, group2tab)

        if (write_out_tables == TRUE) {
            outlieranalysisoutfile = paste(outfilepath, "outlieranalysis_for_",
                groupcombos[groupcombonum,1], "_vs_",
                groupcombos[groupcombonum,2], ".csv", sep="")
            write.table(outtablefin, outlieranalysisoutfile, quote = FALSE,
                        sep = ",", row.names = FALSE)
        }

        outtablelist[[groupcombonum]] = outtablefin
        names(outtablelist)[groupcombonum] = paste("outlieranalysis_for_",
            groupcombos[groupcombonum,1], "_vs_",
            groupcombos[groupcombonum,2], sep="")
    }
    return(outtablelist)
}


## PLOT HEATMAP with metadata, original countdata, and the outlieranalysis to
## pull out significant genes
#' With the grouptablist generated by count_outliers - run through and run a
#'     fisher exact test to get the p.value for the difference in outlier count
#'     for each feature in each of your comparisons
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
#' grouptablist = count_outliers(groupings, outliertab)
#'
#' outlier_analysis_out = outlier_analysis(grouptablist)
#'
#' metatable = sample_annotations
#' counttab = sample_values
#'
#' outlier_heatmap(outlier_analysis_out, analysis_num = NULL,
#'     counttab, metatable, fdrcutoffvalue = 0.1)
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
        GOI = rownames(intable[rowSums(fdrcols < fdrcutoffvalue) >= 1, ])

        if (length(GOI) > 0) {
            print("GOI")
            ## Take the metatable, order it by 1s and then 2s on whatever
            ## comparison we are doing, take comparison columns for plotting
            subsetcounttab = counttab[GOI, rownames(metatable), drop=FALSE]

            annotation1 = annotationlist_builder(metatable)
            outfile1 = paste(outfilepath, "outlieranalysis_heatmap_for_",
                names(outlier_analysis_out)[analysiscount], ".pdf", sep="")
            hm1 = create_heatmap(counttab = subsetcounttab,
                subsetnum = FALSE, colmetatable = metatable,
                colannotationlist = annotation1, rowmetatable = NULL,
                rowannotationlist = NULL, colclusterparam = FALSE,
                rowclusterparam = FALSE, write_out_plot = write_out_plot,
                nameparam = names(outlier_analysis_out)[analysiscount],
                pdfoutfile = outfile1)
            #return(hm1)
            heatmaplist[analysiscount] = list(hm1)
            names(heatmaplist)[analysiscount] = paste0("print_",
                                names(outlier_analysis_out)[analysiscount])
        } else {print(paste("No Significant Outliers for ",
            names(outlier_analysis_out)[analysiscount],
            " at an FDR cut off value of ", fdrcutoffvalue, sep = ""))}
    }
    return(Filter(Negate(is.null), heatmaplist))
}


## COMPLETE BLACKSHEEP FUNCTION
#' Run the entire Blacksheep Function from Start to finish
#' @param counttable the count table to analyze
#' @param metatable the metatable with each column representing a comparison to
#'     perform.
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
#' @param fdrcutoffvalue DEFAULT: 0.1; The FDR value for significance
#' @param outfilepath the full string path to where the file should output to,
#'     DEFAULT is current working directory
#' @return outputs a pdf with the heatmap in the current working directory
#' @keywords outliers
#' @import ComplexHeatmap RColorBrewer circlize
#' @export
#' @examples
#'
#' data("sample_values")
#' data("sample_annotations")
#'
#' blacksheep(counttable = sample_values, metatable = sample_annotations,
#'     analyze_negative_outliers = TRUE, aggregate_features = TRUE,
#'     feature_delineator = "\\.", fdrcutoffvalue = 0.1,
#'     outfilepath = getwd())
blacksheep <- function(counttable, metatable, analyze_negative_outliers = FALSE,
                        aggregate_features = FALSE, feature_delineator = "\\.",
                        fdrcutoffvalue = 0.1, outfilepath = getwd()) {

    ## Use the groupings function to create comparison groups
    groupings = comparison_groupings(metatable)

    ## Create the outlier table and other outlier objects
    reftable_function_out = make_outlier_table(counttable,
                analyze_negative_outliers = analyze_negative_outliers,
                aggregate_features = aggregate_features,
                feature_delineator = "\\.")
    outliertab = reftable_function_out$outliertab
    upperboundtab = reftable_function_out$upperboundtab
    lowerboundtab = reftable_function_out$lowerboundtab
    sampmedtab = reftable_function_out$sampmedtab
    aggposoutlierstab = reftable_function_out$aggposoutlierstab
    aggnegoutlierstab = reftable_function_out$aggnegoutlierstab
    aggposfractiontab = reftable_function_out$aggposfractiontab
    aggnegfractiontab = reftable_function_out$aggnegfractiontab

    if (aggregate_features == FALSE) {
        grouptablist = count_outliers(groupings, outliertab)
        outlier_analysis_out = outlier_analysis(grouptablist)
        pos_outlier_analysis_out = neg_outlier_analysis_out = NULL
        return(list(outlier_analysis = outlier_analysis_out))

        outlier_heatmap(outlier_analysis_out, analysis_num = NULL,
                        counttable, metatable, fdrcutoffvalue = 0.1)

    } else {
        outlier_analysis_out = NULL
        posgrouptablist = count_outliers(groupings, aggposoutlierstab)
        pos_outlier_analysis_out = outlier_analysis(posgrouptablist,
                        outfilepath = paste0(getwd(), "/positive_"))

        outlier_heatmap(outlier_analysis_out = pos_outlier_analysis_out,
                        analysis_num = NULL, counttab = aggposfractiontab,
                        metatable = metatable, fdrcutoffvalue = 0.1,
                        outfilepath = paste0(getwd(), "/positive_"))

        if (analyze_negative_outliers == TRUE) {
            neggrouptablist = count_outliers(groupings, aggnegoutlierstab,
                                            analyze_negative_outliers = TRUE)
            neg_outlier_analysis_out = outlier_analysis(neggrouptablist,
                        outfilepath = paste0(getwd(), "/negative_"))

            outlier_heatmap(neg_outlier_analysis_out, analysis_num = NULL,
                            aggnegfractiontab, metatable, fdrcutoffvalue = 0.1,
                            outfilepath = paste0(getwd(), "/negative_"))
        }

        return(list(positive_outlier_analysis = pos_outlier_analysis_out,
                negative_outlier_analysis = neg_outlier_analysis_out))
    }

}

