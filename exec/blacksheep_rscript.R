!# /usr/bin/rscript

################################
# Name: MacIntosh Cornwell
# Email: mcornwell1957@gmail.com
################################
## Outlier Analysis

## Load in Libraries
print("Beginning Outlier Analysis")
packagelist = c("tools", "argparse", "BlackSheep")
junk <- lapply(packagelist, function(xxx) suppressMessages(
  require(xxx, character.only = TRUE,quietly=TRUE,warn.conflicts = FALSE)))

## ARGPARSER
#####
parse_args <- function() {
  parser <- argparse::ArgumentParser(
    description="BlackSheep Outlier Analysis")
  parser$add_argument('-it', '--intable', help="Path to intable for countdata",
                      type="character", nargs=1)
  parser$add_argument('-im', '--inmetatable', help="Path to intable for metadata",
                      type="character", nargs=1)
  parser$add_argument('-op', '--outpath',
                      help="Path to outpath ex) /path/to/folder/output/",
                      type="character",nargs=1)
  args <- parser$parse_args()
  return (args)
}
args <- parse_args()
if (is.null(args$intable)) {
    stop("Incorrect Input, Enter < Rscript rnaseq_processing.R -h > for help",
         call. = FALSE)}
if (is.null(args$inmetatable)) {
    stop("Incorrect Input, Enter < Rscript rnaseq_processing.R -h > for help",
         call. = FALSE)}
if (is.null(args$outpath)) {
    stop("Incorrect Input, Enter < Rscript rnaseq_processing.R -h > for help",
         call. = FALSE)}

counttabfile = args$intable
metatablefile = args$inmetatable
outfilepath = args$outpath
dir.create(outfilepath, recursive = TRUE, showWarnings = FALSE)

#####


counttab = read.table(normcounttabfile, header = TRUE, row.names = 1,
    sep = ifelse(tools::file_ext(metatablefile)=="txt", "\t", ","),
    check.names = FALSE)
metatable = read.table(metatablefile, header = TRUE, row.names = 1,
    sep = ifelse(tools::file_ext(metatablefile)=="txt", "\t", ","),
    stringsAsFactors = FALSE, check.names = FALSE, na.strings = c(NA, "NA", ""))

## SAFEGUARDS
## First - make sure that count and metatable are working with the same samples
namecheck = intersect(colnames(normcounttab), rownames(metatable))
metatable = metatable[namecheck,]
normcounttab = normcounttab[,namecheck]

## Create groupings based off of metatable
groupings = comparison_groupings(metatable)

## Create outlier table based off of count table
reftable_function_out = make_outlier_table(counttab)
outliertab = reftable_function_out$outliertab

## TABULATE OUTLIER DATA on a per gene basis for each of our subgroups
grouptablist = count_outliers(groupings, outliertab)

## Create a blank list, and fill it with the pvalue for the fisher.test for
## each gene in the upward and downward significance directions
outlier_analysis_out = outlier_analysis(grouptablist,
    analyzenegativeoutliers = FALSE, outfilepath = outfilepath)

## For any significant genes, we then want to add a functionality that will plot
## a heatmap using the metadata and the significant genes
# analysis_num_param = 4
outlier_heatmap(outlier_analysis_out = outlier_analysis_out,
                counttab = counttab,
                metatable = metatable, fdrcutoffvalue = 0.1,
                outfilepath = outfilepath)


