## Internal Functions for Outlier_Analysis function:
split_and_fish <- function(combined_grouptab, sidedness){
    sigval <- ifelse(sidedness=="pos", 1, -1)
    conttab <- matrix(unlist(combined_grouptab), nrow = 2, ncol = 2,
                    dimnames = list(c("0", sigval), c("group1","group2")))
    return(fisher.test(conttab, alternative = "two.sided")$p.value)
}

filter_features <- function(group1tab, group2tab) {
    ## Apply the raw number filter to group1 > 2
    groupprop_filter <- (
        group1tab[,2] / (group1tab[,1] + group1tab[,2])) >
        (group2tab[,2] / (group2tab[,1] + group2tab[,2]))
    groupprop_filter[is.na(groupprop_filter)] <- FALSE
    ## Return a list of genes that pass the raw number filter
    groupprop_filter_features <- rownames(group1tab[groupprop_filter,])
    return(groupprop_filter_features)
}

fraction_filter <- function(fraction_table, groupsamps,
                            fraction_samples_cutoff,
                            groupprop_filter_features) {
    ## apply filter fraction
    groupfractab <- fraction_table[,groupsamps, drop=F]
    groupfractab_select <- groupfractab[rowSums(
        groupfractab!=0, na.rm = TRUE)/ncol(groupfractab) >
            fraction_samples_cutoff,]
    return(rownames(groupfractab_select))
}

fish_to_fdr <- function(fishout, groupprop_filter_features){
    fishgroup <- fishout[rownames(fishout) %in%
                            groupprop_filter_features,, drop = FALSE]
    ## Left in long form to maintain names and labels
    if (nrow(fishgroup) > 1) {fdrvals <- apply(fishgroup, 2, function(x)
        p.adjust(x, method = "BH"))
    } else {
        if(nrow(fishgroup) == 1){
            fdrvals <- fishgroup
        } else {
            fdrvals <- NULL
        }
    }
    return(fdrvals)
}

