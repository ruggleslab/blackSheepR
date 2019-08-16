################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Plotting Functions

## WRITE HEATMAP FUNCTION
#' Plot out a heatmap
#' @param counttab table with counts, samples -x-axis, features -y-axis
#' @param colmetatable the metatable containing information for the columns
#' @param colannotationlist annotation table for columns, based off colmetatable
#' @param colclusterparam cluster the columns?
#' @param rowclusterparam cluster the rows?
#' @param nameparam the title on the heatmap
#' @param write_out_plot DEFAULT: FALSE write out the plot to <pdfoutfile>
#' @param pdfoutfile the full path to the pdf to write out to
#' @return prints a pdf heatmap out to the designated outpath
#' @keywords outliers
#' @import stats grid grDevices ComplexHeatmap circlize
#' @export
#' @examples
#' data("sample_values")
#' counttab = sample_values
#' pdfoutfile = paste(getwd(), "/heatmap_out.pdf", sep = "")
#' nameparam = "testplot"
#'
#' create_heatmap(counttab = counttab,
#'     colmetatable = NULL,
#'     colannotationlist = NULL,colclusterparam = FALSE,
#'     rowclusterparam = FALSE, nameparam, write_out_plot = FALSE,
#'     pdfoutfile= pdfoutfile)
create_heatmap <- function(counttab = counttab,
    colmetatable = NULL, colannotationlist = NULL, colclusterparam = FALSE,
    rowclusterparam = FALSE, nameparam, write_out_plot = FALSE,
    pdfoutfile = pdfoutfile) {

    ## Calc. spearman correlation and use values for column clustering
    if (colclusterparam != FALSE) {
        ## FAILSAFE - Exlcude columns with variance of 0 for plotting
        nearZeroVarcols <- which(vapply(counttab, var) == 0)
        if (length(nearZeroVarcols) > 0) {
            counttab = counttab[,!colnames(counttab) %in%
                                    names(nearZeroVarcols)]
            print(paste("excluding column ", names(nearZeroVarcols),
                " due to zero variance", sep=""))}
        cordata <- cor(counttab, method="spearman")
        coldistance = dist(t(as.matrix(na.omit(cordata))), method = "euclidean")
        colcluster = hclust(coldistance, method = "ward.D2")
        colclusterparam = colcluster
    }

    ## Zscore out counttable
    #maptab = as.matrix(t(apply(counttab, 1, function(x) scale(x))))
    maptab = as.matrix(counttab)

    if (rowclusterparam != FALSE) {
        rowdistance = dist(na.omit(maptab), method = "euclidean")
        rowcluster = hclust(rowdistance, method = "ward.D2")
        rowclusterparam = rowcluster
    }

    ## Build Annotations from our metatable
    hatop = NULL
    if (!is.null(colannotationlist) & !is.null(colmetatable)) {
        temp1 <- vector("list", length(colannotationlist))
        names(temp1) = names(colannotationlist)
        annotlegendlist = lapply(temp1, function(x) x[[1]] = list(title_gp =
            gpar(fontsize = 6, fontface = "bold"), labels_gp=gpar(fontsize=6)))
        ## Param that  keeps a legend if it's annotation is continuous
        ## or has less than 10 discrete terms, otherwise hide the legend
        showlegendparam = unname(unlist(lapply(colannotationlist, function(x) {
            numterms = tryCatch(length(na.omit(unique(x))),
                                error=function(e) NULL)
            is.null(numterms) || numterms <= 10})))
        hatop = HeatmapAnnotation(df = data.frame(colmetatable),
            col = colannotationlist,
            na_col = "white",
            show_annotation_name = TRUE,
            annotation_name_gp = gpar(fontsize = 8, fontface="bold"),
            annotation_name_side = "left",
            simple_anno_size = unit(min(30/length(colannotationlist), 5),"mm"),
            show_legend = showlegendparam,
            annotation_legend_param = annotlegendlist)
    }

    ## Define the Heatmap
    #heatmapcolorparam = colorRamp2(c(-3,0,3), c("blue", "white", "red"))
    lowcolor = midcolor = highcolor = lowvalue = midvalue = highvalue = NULL
    if(max(maptab, na.rm = TRUE) > 0) {
        highcolor = "red"
        highvalue = max(max(maptab, na.rm = TRUE), 1)
    } else {
        highcolor = "white"
        highvalue = 0
    }
    if(min(maptab, na.rm = TRUE) < 0) {
        lowcolor = "blue"
        lowvalue = min(min(maptab, na.rm = TRUE), -1)
    } else {
        lowcolor = "white"
        lowvalue = 0
    }
    if(max(maptab, na.rm = TRUE) > 0 & min(maptab, na.rm = TRUE) < 0) {
        midcolor = "white"
        midvalue = 0
    }
    heatmapcolorparam = colorRamp2(
        c(lowvalue, midvalue, highvalue), c(lowcolor, midcolor, highcolor))
    #heightparam = ifelse(nrow(rnatable) < 20,unit(1,"cm"),NULL)

    ht1 = Heatmap(maptab,
            col = heatmapcolorparam,    ## Define the color scale
            row_title = "Features",                     ## Name the rows
            row_title_gp = gpar(fontsize = 8),
            column_title = nameparam,                   ## Name the columns
            column_title_gp = gpar(fontsize = 8),

            cluster_columns = colclusterparam,          ## Cluster the columns
            cluster_rows = rowclusterparam,             ## Cluster the rows

            show_column_names = ncol(maptab) <=50,      ## Show the Column Names
            column_names_gp = gpar(fontsize = 6),       ## Column Name Size
            show_row_names = nrow(maptab) <=100,        ## Show Row names
            row_names_side = "left",                    ## Place the row names
            row_names_gp = gpar(fontsize=6),

            show_row_dend = nrow(maptab) <=100,         ## Show row dendrogram
            show_column_dend = TRUE,                    ## Show col dendrogram

            heatmap_legend_param = list(title = "Value",
                #legend_height = unit(2, "cm"),
                title_gp = gpar(fontsize = 8, fontface = "bold")),
            top_annotation = hatop,
            height = unit(min((nrow(maptab)/2), 6),"cm"),
            width = unit(min(ncol(maptab), 9),"cm")
    )

    ## Plot out the heatmap
    if (write_out_plot == TRUE) {pdf(file = pdfoutfile, width=11,height=8.5)}
    outplot = draw(ht1, annotation_legend_side = "bottom",
                    padding = unit(c(5, 20, 5, 5), "mm"))
    if (write_out_plot == TRUE) {junk <- dev.off()}
    return(outplot)

}


## Input the metatable and this will build an annotation, if you want to preset
## the colors - then input them in proper list form, either as named list of
## named colors, or with color brewer
#' Create the annotation object for plotting in a heatmap
#' @param metatable the metatable containing information for the columns
#' @param customcolorlist DEFAULT: NULL, enter colorlist to manually set colors
#' @return return the annotation object
#' @keywords outliers
#' @import stats circlize RColorBrewer viridis
#' @export
#' @examples
#' metatable = data.frame(row.names = c("samp1", "samp2", "samp3", "samp4"),
#'     A = c(rep("high", 2), rep("low", 2)), B = seq(1,7,2))
#' customcolorlist = list(A = c("high" = "red", "low" = "blue"),
#'     B = circlize::colorRamp2(seq(-5, 5, length = 3),
#'     RColorBrewer::brewer.pal(3, "Reds")))
#' annotationlist_builder(metatable, customcolorlist)
annotationlist_builder <- function(metatable, customcolorlist = NULL){
    annotlist <- list()
    #colorlist = colors()[grep('gr(a|e)y', colors(), invert = TRUE)]
    ## Create the color list from manual distinct colors, expanded appropriately
    colorlist = rep(c("dodgerblue2", "#E31A1C", "green4","#6A3D9A", "#FF7F00",
                "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
                "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon", "orchid1",
                "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1",
                "yellow4", "yellow3", "darkorange4", "brown"),10)

    colorcount = 0
    for (colnum in seq_len(ncol(metatable))) {

        if (!is.null(customcolorlist) && colnames(metatable[,colnum,drop=FALSE])
            %in% names(customcolorlist)) {
            annotlist[colnames(metatable[,colnum,drop=FALSE])] =
                customcolorlist[colnames(metatable[,colnum,drop=FALSE])]}

        if (!is.null(customcolorlist) &&
            !colnames(metatable[,colnum,drop=FALSE]) %in% names(customcolorlist)
            | is.null(customcolorlist)) {
            if (!is.numeric(metatable[,colnum])) {

                annotlist[[colnum]] = colorlist[(1+colorcount):
                    (colorcount + length(na.omit(unique(metatable[,colnum]))))]
                colorcount = colorcount + length(na.omit(
                    unique(metatable[,colnum])))

                names(annotlist[[colnum]]) = na.omit(unique(metatable[,colnum]))
                names(annotlist)[[colnum]] = colnames(metatable)[colnum]
            }
            else {collist =colorRamp2(seq(min(metatable[,colnum], na.rm = TRUE),
                            max(metatable[,colnum], na.rm = TRUE), length = 3),
                            brewer.pal(3, "Purples"))
            annotlist[[colnum]] = collist
            names(annotlist)[colnum] = colnames(metatable)[colnum]
            }
        }
    }

    return(annotlist)
}

