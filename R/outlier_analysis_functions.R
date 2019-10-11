##########################################
# Name: MacIntosh Cornwell
# Email: macintosh.cornwell@nyulangone.org
##########################################
## blacksheep Script File


## DEVA NORMALIZATION
#' Normalization of data to prepare for deva. Uses a Median of Ratio method
#' followed by a log2 transformation.
#'
#' @param intable table with samples along the columns and features along the
#'     rows.
#' @param method DEFAULT: "MoR-log"; Method by which to normalize data in
#'     preparation for deva. Options are <"MoR-log", "MoR", "log">. Where "MoR"
#'     refers to the Median of ratio's. The "log" transformation is necessary
#'     to compress heavily skewed data and allow for proper detection.
#'     "MoR-log" as the default will perform MoR followed by a log2 transform.
#' @return A normalized table for input into deva
#' @keywords outliers blacksheepr deva
#' @import stats pasilla
#' @export
#' @examples
#' library(pasilla)
#' pasCts <- system.file("extdata",
#'     "pasilla_gene_counts.tsv", package="pasilla")
#' cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
#' norm_cts <- deva_normalization(cts, method = "MoR-log")
deva_normalization <- function(intable, method = "MoR-log"){
    ## Create empty outtablelist
    normmethods <- unlist(strsplit(method, "-"))

    ## Perform MoR to data if called for
    if ("MoR" %in% normmethods) {
        pseudoref <- log(exp(rowMeans(log(intable))))
        norm_factor <- apply(intable, 2, function(x) {
            exp(median((log(x) - pseudoref)[is.finite(pseudoref) & x > 0]))})
        norm_factor <- norm_factor/exp(mean(log(norm_factor)))
        outtable <- as.matrix(intable) %*% diag(1 / norm_factor)
        colnames(outtable) <- colnames(intable)
        if ("log" %in% normmethods) {
            outtable <- log2(outtable)
        }
    ## If no MoR, see if we still want to log transform
    } else {
        if ("log" %in% normmethods) {
            outtable <- log2(intable)
        }
    }
    return(outtable)
}


## MAKE COMPARISON COLUMNS
#' Utility function that will take in columns with several subcategories,
#' and output several columns each with binary classifications.
#' ex) col1: A,B,C >> colA: A,notA,notA; colB: notB,B,notB; colC: notC,notC,C
#'
#' @param intable table where each column has more than one subcategory, can
#'     be multiple columns
#' @return an expanded table with each of the columns as a binary labeling of
#'     each subcategory.
#' @keywords outliers blacksheepr deva
#' @import stats
#' @export
#' @examples
#' data("sample_annotationdata")
#' new_comparisons <- make_comparison_columns(
#'     sample_annotationdata[,1,drop=FALSE])
make_comparison_columns <- function(intable){
    ## Create empty outtablelist
    outtablelist <- list()
    ## Run through each inputted column
    for (i in seq_len(ncol(intable))) {
        ## Pull out the rownames and column names to save for later
        featurenames <- rownames(intable)
        categoryname <- colnames(intable)[i]

        ## Select the column we are working with
        subtable <- intable[,i, drop=FALSE]
        subtable <- apply(subtable, 2, as.character)

        ## Create duplicate columns - one for each subcategory
        new_comp_tab <- matrix(rep(subtable, length(unique(subtable[,1]))),
                    nrow = nrow(intable), ncol = length(unique(subtable[,1])),
                    dimnames = list(featurenames, unique(subtable[,1])))

        ## For each subcategory, declare it as binary
        for ( ii in seq_len(ncol(new_comp_tab))) {
            new_comp_tab[new_comp_tab[,ii, drop=FALSE] !=
                            colnames(new_comp_tab[,ii, drop=FALSE]),ii] <-
                paste0("not_", colnames(new_comp_tab[,ii,drop=FALSE]))
        }

        ## Add back the column name, and then save out
        colnames(new_comp_tab) <- paste0(categoryname, "_",
                                        colnames(new_comp_tab))
        outtablelist[[i]] <- new_comp_tab

    }
    output <- do.call(cbind, outtablelist)
    return(output)
}


## EXTRACT COMPARISONS
#' Create all of the groups based on the input metadata
#'
#' @param comptable table where each column will have comparisons drawn from it
#' @return a list with each of the groups as an entry in the list
#'     NOTE - this list will be ncol*2 long where ncol is the number comparisons
#' @keywords outliers blacksheepr deva
#' @import stats
#' @export
#' @examples
#' data("sample_annotationdata")
#' groupings <- comparison_groupings(sample_annotationdata)
comparison_groupings <- function(comptable) {
    ## Initate blank list
    groupings <- list()

    ## FAILSAFE - to insure that these are charactors not factors
    comptable <- data.frame(comptable)
    samplenames <- rownames(comptable)
    comptable <- data.frame(lapply(comptable,as.character),
                                stringsAsFactors = FALSE)
    rownames(comptable) <- samplenames
    for (i in seq_len(ncol(comptable))) {
        # Select column and create comparison based off entries in column
        subsamp <- comptable[,i, drop=FALSE]
        subsamp[subsamp == ""] <- NA ## Failsafe to turn blank cells into NA
        subcats <- unique(na.omit(subsamp[,1]))

        ## Fill the i'th entry in list with samples that fall into subcategory
        groupings[[i]] <- rownames(na.omit(subsamp[subsamp == subcats[1],,
                                                drop=FALSE]))
        names(groupings)[i] <- paste0(colnames(subsamp), "__", subcats[1])

        ## fill the i+num comps entry with samples in that subcategory
        groupings[[i+ncol(comptable)]] <- rownames(na.omit(
            subsamp[subsamp == subcats[2],,drop=FALSE]))
        names(groupings)[i+ncol(comptable)] <- paste0(
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
#'
#' @param intable table with all of the inputted information, samples along the
#'     x-axis, features along the y-axis
#' @param analyze_negative_outliers DEFAULT: FALSE; Toggle the analysis of
#'     outliers in the negative direction. Will lead to the output of
#'     the outlier table containing "-1" values, in addition to negative outputs
#'     for boundaries and aggregate tables (if applicable)
#' @return a list with varied sections depending on parameters:
#'     $outliertab - table converted to outlier form with 0s, 1s, and -1s,
#'     $upperboundtab - list of upper boundaries for outliers
#'     $lowerboundtab - list of lower boundaries of outliers
#'     $sampmedtab - list of median value per feature
#' @keywords outliers blacksheepr deva
#' @export
#' @examples
#' data("sample_phosphodata")
#' reftable_function_out <- make_outlier_table(sample_phosphodata[1:1000,],
#'     analyze_negative_outliers = FALSE)
#' outliertab <- reftable_function_out$outliertab
#' upperboundtab <- reftable_function_out$upperboundtab
#' lowerboundtab <- reftable_function_out$lowerboundtab
#' sampmedtab <- reftable_function_out$sampmedtab
make_outlier_table <- function(intable, analyze_negative_outliers = FALSE){

    ## Define blank out lists
    outlist <- uppboundlist <- lowboundlist <- sampmedlist <- list()
    for (intablegenecount in seq_len(nrow(intable))) {

        ## Save the sampname, then subset the table to just the one gene
        genename <- rownames(intable)[intablegenecount]
        sampdat <- as.data.frame(t(intable[intablegenecount,, drop=FALSE]))
        sampdat$categ <- 0
        sampdat[is.na(sampdat[,1]),2] <- NA
        colnames(sampdat) <- c(genename, genename)

        ## Set upper and lower bounds based off of the IQR and med
        sampmedlist[[intablegenecount]] <- sampmed <-
            median(sampdat[,1],na.rm=TRUE)
        names(sampmedlist)[intablegenecount] <- genename

        if (analyze_negative_outliers != TRUE) {
            uppboundlist[[intablegenecount]] <- uppbound <-
                sampmed + 1.5*IQR(sampdat[,1], na.rm = TRUE)
            names(uppboundlist)[intablegenecount] <- genename
            sampdat[!is.na(sampdat[,1]) & sampdat[,1] > uppbound,2] <- 1
            lowboundlist <- NULL
        } else {
            lowboundlist[[intablegenecount]] <- lowbound <-
                sampmed - 1.5*IQR(sampdat[,1], na.rm = TRUE)
            names(lowboundlist)[intablegenecount] <- genename
            sampdat[!is.na(sampdat[,1]) & sampdat[,1] < lowbound,2] <- -1
            uppboundlist <- NULL
        }

        ## Save what we just made - which is n by 2 table with the gene in one
        ## col and the status of the sample in the other
        outlist[[intablegenecount]] <- sampdat[,2, drop=FALSE]

    }
    ## Combine all of our lists
    outliertab = as.data.frame(t(do.call(cbind, outlist)))
    if (analyze_negative_outliers == TRUE) {
        negativeoutlierlist <- list(lowerboundtab =
                                        do.call(rbind, lowboundlist))
        positiveoutlierlist <- NULL
    } else {
        positiveoutlierlist <- list(upperboundtab =
                                        do.call(rbind, uppboundlist))
        negativeoutlierlist <- NULL
    }
    sampmedtab = do.call(rbind, sampmedlist)

    output = c(list(outliertab = outliertab, sampmedtab = sampmedtab),
                positiveoutlierlist, negativeoutlierlist)

    ## Return the outputted values
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
#' @keywords outliers blacksheepr deva
#' @export
#' @examples
#'
#' data("sample_phosphodata")
#' reftable_function_out <- make_outlier_table(sample_phosphodata[1:1000,])
#' outliertab <- reftable_function_out$outliertab
#'
#' data("sample_annotationdata")
#' groupings <- comparison_groupings(sample_annotationdata)
#'
#' count_outliers_out <- count_outliers(groupings, outliertab,
#'     aggregate_features = FALSE)
#' grouptablist <- count_outliers_out$grouptablist
#' fractiontab <- count_outliers_out$fractiontab
count_outliers <- function(groupings, outliertab,
    aggregate_features = FALSE, feature_delineator = "\\.") {

    ## Define empty starting list
    grouptablist <- list()
    ## Set factor depending on analysis - pos/neg - detected automatically
    outliervalue <- unique(unlist(outliertab))[unique(unlist(outliertab)) != 0 &
                                    !is.na(unique(unlist(outliertab)))]

    ## If aggregate is true, split on delineator, and output an outlier and
    ## nonoutlier aggregate table
    if (aggregate_features == TRUE) {
        feature_labels <- do.call(rbind, strsplit(sub(feature_delineator,
                                "xyz", rownames(outliertab)), split = "xyz"))
        colnames(feature_labels) <- c("primary_feature", "secondary_feature")
        aggtabin <- cbind(feature_labels, outliertab)

        aggoutliertab <- aggregate(aggtabin[,c(3:ncol(aggtabin))],
                    by = list(primary_feature = aggtabin[,"primary_feature"]),
                    function(x) sum(x==outliervalue, na.rm = TRUE))
        aggoutliertab <- data.frame(aggoutliertab[,-1],
                        row.names = aggoutliertab[,1], check.names = FALSE)

        aggnonoutliertab <- aggregate(aggtabin[,c(3:ncol(aggtabin))],
                    by = list(primary_feature = aggtabin[,"primary_feature"]),
                    function(x) sum(x!=outliervalue, na.rm = TRUE))
        aggnonoutliertab <- data.frame(aggnonoutliertab[,-1],
                        row.names = aggnonoutliertab[,1], check.names = FALSE)

        ## Create fraction table in similar manner
        fractiontab <- aggregate(aggtabin[,c(3:ncol(aggtabin))],
            by = list(primary_feature = aggtabin[,"primary_feature"]),
                function(x) sum(x[x==outliervalue],na.rm = TRUE)/sum(!is.na(x)))
        fractiontab <- data.frame(fractiontab[,-1],
                        row.names=fractiontab[,1], check.names = FALSE)
    } else {
        ## If no aggregation, then the fractiontab is the same as the outliertab
        if (outliervalue == 1) {
            fractiontab <- outliertab
        } else {
            fractiontab <- -outliertab
        }
    }

    ## Perform counting in desired direction
    for (groupcount in seq_len(length(groupings))) {
        subgroup <- outliertab[,colnames(outliertab) %in%
                                        groupings[[groupcount]], drop=FALSE]
        if (aggregate_features == FALSE) {
            grouptablist[[groupcount]] <- list()
            grouptablist[[groupcount]][[1]] <- t(apply(subgroup,MARGIN=1,
                        function(x) table(factor(x, levels=c(0,outliervalue)))))
            grouptablist[[groupcount]][[2]] <- groupings[[groupcount]]
            names(grouptablist[[groupcount]]) <- c("feature_counts", "samples")

        } else {
            aggoutliercount <- cbind.data.frame(
                primary_feature = rownames(aggoutliertab),
                outlier_count = rowSums(
                                aggoutliertab[,groupings[[groupcount]]]))

            aggnonoutliercount <- cbind.data.frame(
                primary_feature = rownames(aggnonoutliertab),
                outlier_count = rowSums(
                                aggnonoutliertab[,groupings[[groupcount]]]))

            aggcounttab <- merge(aggnonoutliercount, aggoutliercount,
                                by = "primary_feature")
            rownames(aggcounttab) <- aggcounttab[,1]
            colnames(aggcounttab) <- c("primary_feature", 0, outliervalue)

            ## Save out the table to the first entry in the sublist
            grouptablist[[groupcount]] <- list()
            grouptablist[[groupcount]][[1]] <- aggcounttab[,c(2,3)]

            ## Also saving the samples as an additional entry in the sublist
            grouptablist[[groupcount]][[2]] <- groupings[[groupcount]]
            names(grouptablist[[groupcount]]) <- c("feature_counts", "samples")

        }
        names(grouptablist)[groupcount] <- names(groupings)[groupcount]
    }
    ## Define our output with different lists depending on input params
    if (aggregate_features == TRUE){
        aggregateoutlist <- list(aggoutliertab = aggoutliertab)
    } else {
        aggregateoutlist <- NULL
    }
    output <- c(list(grouptablist = grouptablist, fractiontab = fractiontab),
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
#'     write_out_tables = FALSE, outfilepath = tempdir())
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
#'     DEFAULT is a tempdir()
#' @return the analysis table with p.value, fdr, and raw data per comparison
#' @keywords outliers blacksheepr deva
#' @import stats utils
#' @export
#' @examples
#'
#' data("sample_phosphodata")
#' reftable_function_out <- make_outlier_table(sample_phosphodata[1:1000,])
#' outliertab <- reftable_function_out$outliertab
#'
#' data("sample_annotationdata")
#' groupings <- comparison_groupings(sample_annotationdata)
#'
#' count_outliers_out <- count_outliers(groupings, outliertab,
#'     aggregate_features = FALSE)
#' grouptablist <- count_outliers_out$grouptablist
#' fractiontab <- count_outliers_out$fractiontab
#'
#' outlier_analysis_out <- outlier_analysis(grouptablist,
#'     fraction_table = fractiontab)
outlier_analysis <- function(grouptablist,
                        fraction_table = NULL, fraction_samples_cutoff = 0.3,
                        write_out_tables = FALSE, outfilepath = tempdir()) {
    ## Define blank starting lists
    outtablelist <- list()
    ## Create comparison matrix
    groupcombos <- matrix(names(grouptablist),
                        nrow=(length(grouptablist)/2), ncol = 2)
    for (groupcombonum in seq_len(nrow(groupcombos))) {
        ## Pull out the two groups of interest for the run
        group1label <- groupcombos[groupcombonum,1]
        group2label <- groupcombos[groupcombonum,2]
        group1tab <- as.data.frame(grouptablist[[group1label]][[1]])
        group2tab <- as.data.frame(grouptablist[[group2label]][[1]])
        group1samps <- grouptablist[[group1label]][[2]]
        group2samps <- grouptablist[[group2label]][[2]]

        ## Vectorized split_and_fish - extracting pval for each comparison
        comptabtemp <- merge(group2tab, group1tab, by = "row.names")
        comptab <- data.frame(comptabtemp, row.names = comptabtemp[,1])[,2:5]

        sidedparam <- ifelse("1" %in% colnames(group1tab), "pos", "neg")
        statgroup <- apply(comptab, 1, split_and_fish, sidedness = sidedparam)

        fishout <- cbind(statgroup, statgroup)
        colnames(fishout)[c(1,2)] <- paste0("pval_more_", sidedparam,
                            "_outliers_in_", groupcombos[groupcombonum,c(1,2)])

        #### RAW NUMBER FILTER
        ## Filter to only select features that already have a proportion of
        ## outliers greater in the ingroup
        group1prop_filter_features <- filter_features(group1tab, group2tab)
        group2prop_filter_features <- filter_features(group2tab, group1tab)
        prop_filter_features <- union(group1prop_filter_features,
                                    group2prop_filter_features)

        ## If we are applying the fraction table, then apply this filter too
        if (!is.null(fraction_table)){
            group1fractab_features <- fraction_filter(fraction_table,
                                group1samps, fraction_samples_cutoff,
                                group1prop_filter_features)

            group2fractab_features <- fraction_filter(fraction_table,
                                group2samps, fraction_samples_cutoff,
                                group2prop_filter_features)

            fraction_selected_genes <- union(group1fractab_features,
                                            group2fractab_features)

        } else {fraction_selected_genes <- rownames(fishout)}

        fishout <- fishout[rownames(fishout) %in% intersect(
            fraction_selected_genes, prop_filter_features),,drop=FALSE]

        ## Apply the fdr in both directions separately
        fdrvals1 <- fish_to_fdr(fishout, group1prop_filter_features)
        fdrvals2 <- fish_to_fdr(fishout, group2prop_filter_features)
        fdrvals <- rbind(fdrvals1, fdrvals2)
        if (is.null(fdrvals)){next} # break to skip iteration if no sig genes
        colnames(fdrvals) <- gsub(pattern = "pval", replacement = "fdr",
                                colnames(fdrvals))

        ## SAVE OUT DATA
        colnames(group1tab) <- paste0(group1label, "_", paste0(c("non", ""),
                                                    sidedparam, "outlier"))
        colnames(group2tab) <- paste0(group2label, "_", paste0(c("non", ""),
                                                    sidedparam, "outlier"))

        ## Write out the final table
        outtablefin <- Reduce(function(dtf1, dtf2)
            merge(dtf1, dtf2,by = "gene"), lapply(list(fishout, fdrvals,
                group1tab[fraction_selected_genes,],
                group2tab[fraction_selected_genes,]),
                function(x) data.frame(x, gene = row.names(x))))

        ## Format the final table to only have the pvalues for the analysis that
        ## leads to more in ingroup/outgroup
        outtablefin[outtablefin[,7]/(outtablefin[,6] + outtablefin[,7]) >
            outtablefin[,9]/(outtablefin[,8] + outtablefin[,9]),c(3,5)] <- ""
        outtablefin[outtablefin[,7]/(outtablefin[,6] + outtablefin[,7]) <
            outtablefin[,9]/(outtablefin[,8] + outtablefin[,9]),c(2,4)] <- ""

        ## Write out the tables
        if (write_out_tables == TRUE) {
            outlieranalysisoutfile <- paste(outfilepath, "outlieranalysis_for_",
                group1label, "_vs_", group2label, ".csv", sep="")
            write.table(outtablefin, outlieranalysisoutfile, quote = FALSE,
                        sep = ",", row.names = FALSE)
        }

        ## Save out results in a list
        outtablelist[[groupcombonum]] <- outtablefin
        names(outtablelist)[groupcombonum] = paste("outlieranalysis_for_",
            group1label, "_vs_", group2label, sep="")
    }
    return(Filter(Negate(is.null), outtablelist))
}


## PLOT HEATMAP with metadata, original countdata, and the outlieranalysis to
## pull out significant genes
#' With the grouptablist generated by count_outliers - run through and run a
#'     fisher exact test to get the p.value for the difference in outlier count
#'     for each feature in each of your comparisons
#' @usage outlier_heatmap(outlier_analysis_out, analysis_num = NULL,
#'     counttab, metatable, fdrcutoffvalue = 0.1)
#' @param outlier_analysis_out the full outlier_analysis data objet
#' @param analysis_num DEFAULT: NULL; if you only want to plot the heatmap for
#'     a particular analysis, enter number of that analysis
#' @param counttab the raw data before outlier analysis
#' @param metatable the complete metatable that was used to generate the
#'     comparisons, will be used for annotation of the heatmap
#' @param fdrcutoffvalue DEFAULT: 0.1; The FDR value for significance
#' @return outputs a pdf with the heatmap in the current working directory
#' @keywords outliers blacksheepr deva
#' @import ComplexHeatmap RColorBrewer circlize
#' @export
#' @examples
#'
#' data("sample_phosphodata")
#' reftable_function_out <- make_outlier_table(sample_phosphodata[1:1000,])
#' outliertab <- reftable_function_out$outliertab
#'
#' data("sample_annotationdata")
#' groupings <- comparison_groupings(sample_annotationdata)
#'
#' count_outliers_out <- count_outliers(groupings, outliertab,
#'     aggregate_features = FALSE)
#' grouptablist <- count_outliers_out$grouptablist
#' fractiontab <- count_outliers_out$fractiontab
#'
#' outlier_analysis_out <- outlier_analysis(grouptablist,
#'     fraction_table = fractiontab)
#'
#' metatable <- sample_annotationdata
#' counttab <- sample_phosphodata
#'
#' hm1 <- outlier_heatmap(outlier_analysis_out, analysis_num = NULL,
#'     fractiontab, metatable, fdrcutoffvalue = 0.1)
outlier_heatmap <- function(outlier_analysis_out, analysis_num = NULL, counttab,
                            metatable, fdrcutoffvalue = 0.1) {
    ## Define the start and stops of the foorloop - if they put in a single
    ## analysis to do - it will only output that, otherwise, will loop over all
    startcount <- ifelse(is.null(analysis_num), 1, analysis_num)
    endcount <- ifelse(is.null(analysis_num), length(outlier_analysis_out),
                    analysis_num)
    heatmaplist <- NULL
    for (analysiscount in startcount:endcount) {
        intable <- outlier_analysis_out[[analysiscount]]

        ## Select the column from the outlier_analysis_out that contain "fdr"
        ## and then grab rownames for columns that have sig value
        fdrcols <- intable[ ,grepl( "fdr", colnames(intable))]
        GOI <- intable[rowSums(fdrcols < fdrcutoffvalue) == 2,1]

        if (length(GOI) > 0) {
            ## Take the metatable, order it by 1s and then 2s on whatever
            ## comparison we are doing, take comparison columns for plotting
            ## Added in as.character as a failsafe
            subsetcounttab <- counttab[as.character(GOI),
                                        rownames(metatable), drop=FALSE]

            annotation1 <- annotationlist_builder(metatable)
            hm1 <- create_heatmap(counttab = subsetcounttab,
                colmetatable = metatable, colannotationlist = annotation1,
                colclusterparam = FALSE, rowclusterparam = FALSE,
                nameparam = names(outlier_analysis_out)[analysiscount])
            #return(hm1)
            heatmaplist[analysiscount] <- list(hm1)
            names(heatmaplist)[analysiscount] <- paste0("print_",
                                names(outlier_analysis_out)[analysiscount])
        }
    }
    return(Filter(Negate(is.null), heatmaplist))
}


## COMPLETE BLACKSHEEP FUNCTION
#' Run the entire blacksheep Function from Start to finish
#' @usage deva(se, analyze_negative_outliers = FALSE,
#'     aggregate_features = FALSE, feature_delineator = "\\\\.",
#'     fraction_samples_cutoff = 0.3, fdrcutoffvalue = 0.1)
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
#' @return outputs the full output of deva, including the analysis tables, the
#'     heatmaps for the analyses, the fraction table showing the fraction of
#'     outliers per sample, and the median and boundary values that together
#'     comprise the outlier boundary
#' @keywords outliers blacksheepr deva
#' @import ComplexHeatmap RColorBrewer circlize
#' @rawNamespace import(SummarizedExperiment, except = c(start, end))
#' @export
#' @examples
#'
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#' data("sample_phosphodata")
#' data("sample_annotationdata")
#'
#' se <- SummarizedExperiment(
#'     assays = list(counts = as.matrix(sample_phosphodata[1:1000,])),
#'     colData = DataFrame(sample_annotationdata))
#'
#' deva(se = se,
#'     analyze_negative_outliers = FALSE, aggregate_features = FALSE,
#'     feature_delineator = "-", fraction_samples_cutoff = 0.3,
#'     fdrcutoffvalue = 0.1)
#'
#'
deva <- function(se, analyze_negative_outliers = FALSE,
                        aggregate_features = FALSE, feature_delineator = "\\.",
                        fraction_samples_cutoff = 0.3, fdrcutoffvalue = 0.1) {

    ## Extracting information from the SummarizedExperiment
    counttable <- assays(se)[[1]]
    metatable <- colData(se)

    ## Use the groupings function to create comparison groups
    groupings <- comparison_groupings(metatable)

    ## Create the outlier table and other outlier objects
    reftable_function_out <- make_outlier_table(counttable,
                analyze_negative_outliers = analyze_negative_outliers)
    outliertab <- reftable_function_out$outliertab
    upperboundtab <- reftable_function_out$upperboundtab
    lowerboundtab <- reftable_function_out$lowerboundtab
    sampmedtab <- reftable_function_out$sampmedtab

    ## Count the number of outliers
    count_outliers_out <- count_outliers(groupings, outliertab,
                                aggregate_features = aggregate_features,
                                feature_delineator = feature_delineator)
    grouptablist <- count_outliers_out$grouptablist
    fractiontab <- count_outliers_out$fractiontab

    ## Perform outlier analysis
    outlier_analysis_out <- outlier_analysis(
        grouptablist = grouptablist,
        fraction_table = fractiontab,
        fraction_samples_cutoff = fraction_samples_cutoff)

    ## Expanded out Code to for loop to sort the metatable for each comp
    hm1list <- list()
    for (i in seq_len(length(outlier_analysis_out))){
        plottable <- metatable[do.call(order, c(decreasing = TRUE,
                            data.frame(metatable[,c(i,
                                setdiff(seq_len(ncol(metatable)),i))]))),]
        hm1list[[i]] <- outlier_heatmap(
            outlier_analysis_out = outlier_analysis_out, analysis_num = i,
            counttab = fractiontab, metatable = plottable, fdrcutoffvalue)
    }
    hm1 <- unlist(hm1list)

    return(list(outlier_analysis_out = outlier_analysis_out,
                significant_heatmaps = hm1, fraction_table = fractiontab,
                median_value_table = sampmedtab,
                outlier_boundary_table = c(upperboundtab, lowerboundtab)))
}

## HELPER FUNCTIONS:

## Deva Results Summary
#' Utility function that allows easier grabbing of data
#'
#' @param deva_out output from the deva function
#' @param ID The keyword to search through analyses and grab desired output
#' @param type <"table" | "heatmap" | "fraction_table" | "median" | "boundary">
#'     to return the desirted analysis type
#' @return desired subset of analysis from deva
#' @keywords outliers blacksheepr deva
#' @export
#' @examples
#'
#' suppressPackageStartupMessages(library(SummarizedExperiment))
#' data("sample_phosphodata")
#' data("sample_annotationdata")
#'
#' se = SummarizedExperiment(
#'     assays = list(counts = as.matrix(sample_phosphodata[1:1000,])),
#'     colData = DataFrame(sample_annotationdata))
#'
#' deva_out = deva(se = se,
#'     analyze_negative_outliers = FALSE, aggregate_features = TRUE,
#'     feature_delineator = "-", fraction_samples_cutoff = 0.3,
#'     fdrcutoffvalue = 0.1)
#'
#' deva_results(deva_out, ID = "outlieranalysis", type = "table")

deva_results <- function(deva_out, ID = NULL, type = NULL) {
    # If just deva_out is input - return names of analyses performed
    if ((is.null(ID))) {
        if(is.null(type)){
            out1 <- names(deva_out)
        } else {
            if(type == "table"){out1 <- names(deva_out[[1]])}
            if(type == "heatmap"){out1 <- names(deva_out[[2]])}
            if(type == "fraction_table"){out1 <- deva_out[[3]]}
            if(type == "median"){out1 <- deva_out[[4]]}
            if(type == "boundary"){out1 <- deva_out[[5]]}
        }
    } else {
        ## Determine what type of object to output
        if(is.null(type)) {stop("Please input data type to output
            <\"table\" | \"heatmap\" | \"fraction_table\" |
            \"median\" | \"boundary\">")}

        if(type == "table"){datacat <- 1}
        if(type == "heatmap"){datacat <- 2}

        ## If the table of heatmap - named analyses, then search for name
        if(datacat == 1|2) {
            grab <- grep(pattern = ID, names(deva_out[[datacat]]))

            if (length(grab) == 1) {out1 <- deva_out[[datacat]][[grab]]}
            if (length(grab) == 0) {stop("No analyses have name with \"", ID,
                    "\" Please check output names <deva_results(deva_out)>")}
            if (length(grab) > 1) {
                warning("Multiple analyses matched ", ID,
                                ", outputting named list of selected analyses")
                out1 <- deva_out[[datacat]][grab]
        }} else {
        ## If just looking for the fraction_table, median, or boundary table
            out1 <- deva_out[[datacat]]
        }

    }

    return(out1)
}


