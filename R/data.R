#' sample_annotationdata
#'
#' Example annotation data for Outlier analysis.
#' This example data is a subset of the data used in the CPTAC3 Breast Cancer
#' exploration study: (doi: 10.1038/nature18003).
#' Each row corresponds to a sample and each column is an binary annotation for
#' that sample.
#'
#' @format A data frame with 76 rows and 6 variables:
#' \describe{
#'     \item{PAM50_Her2}{The binary PAM50 Her2 classification for each sample}
#'     \item{PAM50_Basal}{The binary PAM50 Basal classification for each sample}
#'     \item{PAM50_LumA}{The binary PAM50 LumA classification for each sample}
#'     \item{PAM50_LumB}{The binary PAM50 LumB classification for each sample}
#'     \item{ER_Status}{The ER Status classification for each sample}
#'     \item{PR_Status}{The PR Status classification for each sample}
#'     ...
#' }
#' @source \url{https://cptac-data-portal.georgetown.edu/cptac/s/S029}
"sample_annotationdata"


#' sample_phosphodata
#'
#' Example phosphoprotein data for Outlier analysis
#' This example data is a subset of the data used in the CPTAC3 Breast Cancer
#' exploration study: (doi: 10.1038/nature18003).
#' Each row corresponds to a phosphoprotein site, and each column is a sample.
#' The values within the table are normalized massspec phosphoprotein values.
#'
#' @format A data frame with 15532 rows and 76 variables:
#' \describe{
#'     \item{TCGA-A2-A0CM}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0D2}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0EQ}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0EV}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0EX}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0EY}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0SW}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0SX}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0T3}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0T6}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0YC}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0YD}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0YF}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0YG}{phosphoprotein levels for each gene}
#'     \item{TCGA-A2-A0YM}{phosphoprotein levels for each gene}
#'     \item{TCGA-A7-A0CE}{phosphoprotein levels for each gene}
#'     \item{TCGA-A7-A0CJ}{phosphoprotein levels for each gene}
#'     \item{TCGA-A7-A13F}{phosphoprotein levels for each gene}
#'     \item{TCGA-A8-A06N}{phosphoprotein levels for each gene}
#'     \item{TCGA-A8-A06Z}{phosphoprotein levels for each gene}
#'     \item{TCGA-A8-A076}{phosphoprotein levels for each gene}
#'     \item{TCGA-A8-A079}{phosphoprotein levels for each gene}
#'     \item{TCGA-A8-A08Z}{phosphoprotein levels for each gene}
#'     \item{TCGA-A8-A09G}{phosphoprotein levels for each gene}
#'     \item{TCGA-AN-A04A}{phosphoprotein levels for each gene}
#'     \item{TCGA-AN-A0AJ}{phosphoprotein levels for each gene}
#'     \item{TCGA-AN-A0AL}{phosphoprotein levels for each gene}
#'     \item{TCGA-AN-A0AM}{phosphoprotein levels for each gene}
#'     \item{TCGA-AN-A0FK}{phosphoprotein levels for each gene}
#'     \item{TCGA-AN-A0FL}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A03O}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A0J6}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A0J9}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A0JC}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A0JE}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A0JJ}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A0JL}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A0JM}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A126}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A12B}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A12D}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A12E}{phosphoprotein levels for each gene}
#'     \item{TCGA-AO-A12F}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A0TR}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A0TT}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A0TV}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A0TX}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A0U4}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A1AP}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A1AS}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A1AV}{phosphoprotein levels for each gene}
#'     \item{TCGA-AR-A1AW}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A0AV}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A0BV}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A0C1}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A0C7}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A0DD}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A0DG}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A0E1}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A0E9}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A18N}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A18Q}{phosphoprotein levels for each gene}
#'     \item{TCGA-BH-A18U}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A12L}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A12T}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A12U}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A12V}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A12Z}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A130}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A131}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A134}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A135}{phosphoprotein levels for each gene}
#'     \item{TCGA-C8-A138}{phosphoprotein levels for each gene}
#'     \item{TCGA-D8-A142}{phosphoprotein levels for each gene}
#'     \item{TCGA-E2-A154}{phosphoprotein levels for each gene}
#'     \item{TCGA-E2-A158}{phosphoprotein levels for each gene}
#' }
#' @source \url{https://cptac-data-portal.georgetown.edu/cptac/s/S029}
"sample_phosphodata"


#' sample_rnadata
#'
#' Example RNA data for Outlier analysis
#' This example data is a subset of the data used in the CPTAC3 Breast Cancer
#' exploration study: (doi: 10.1038/nature18003).
#' Each row corresponds to a gene, and each column is a sample. The values
#' within the table are normalized transcript counts.
#'
#' @format A data frame with 4317 rows and 76 variables:
#' \describe{
#'     \item{TCGA-A2-A0CM}{RNA levels for each gene}
#'     \item{TCGA-A2-A0D2}{RNA levels for each gene}
#'     \item{TCGA-A2-A0EQ}{RNA levels for each gene}
#'     \item{TCGA-A2-A0EV}{RNA levels for each gene}
#'     \item{TCGA-A2-A0EX}{RNA levels for each gene}
#'     \item{TCGA-A2-A0EY}{RNA levels for each gene}
#'     \item{TCGA-A2-A0SW}{RNA levels for each gene}
#'     \item{TCGA-A2-A0SX}{RNA levels for each gene}
#'     \item{TCGA-A2-A0T3}{RNA levels for each gene}
#'     \item{TCGA-A2-A0T6}{RNA levels for each gene}
#'     \item{TCGA-A2-A0YC}{RNA levels for each gene}
#'     \item{TCGA-A2-A0YD}{RNA levels for each gene}
#'     \item{TCGA-A2-A0YF}{RNA levels for each gene}
#'     \item{TCGA-A2-A0YG}{RNA levels for each gene}
#'     \item{TCGA-A2-A0YM}{RNA levels for each gene}
#'     \item{TCGA-A7-A0CE}{RNA levels for each gene}
#'     \item{TCGA-A7-A0CJ}{RNA levels for each gene}
#'     \item{TCGA-A7-A13F}{RNA levels for each gene}
#'     \item{TCGA-A8-A06N}{RNA levels for each gene}
#'     \item{TCGA-A8-A06Z}{RNA levels for each gene}
#'     \item{TCGA-A8-A076}{RNA levels for each gene}
#'     \item{TCGA-A8-A079}{RNA levels for each gene}
#'     \item{TCGA-A8-A08Z}{RNA levels for each gene}
#'     \item{TCGA-A8-A09G}{RNA levels for each gene}
#'     \item{TCGA-AN-A04A}{RNA levels for each gene}
#'     \item{TCGA-AN-A0AJ}{RNA levels for each gene}
#'     \item{TCGA-AN-A0AL}{RNA levels for each gene}
#'     \item{TCGA-AN-A0AM}{RNA levels for each gene}
#'     \item{TCGA-AN-A0FK}{RNA levels for each gene}
#'     \item{TCGA-AN-A0FL}{RNA levels for each gene}
#'     \item{TCGA-AO-A03O}{RNA levels for each gene}
#'     \item{TCGA-AO-A0J6}{RNA levels for each gene}
#'     \item{TCGA-AO-A0J9}{RNA levels for each gene}
#'     \item{TCGA-AO-A0JC}{RNA levels for each gene}
#'     \item{TCGA-AO-A0JE}{RNA levels for each gene}
#'     \item{TCGA-AO-A0JJ}{RNA levels for each gene}
#'     \item{TCGA-AO-A0JL}{RNA levels for each gene}
#'     \item{TCGA-AO-A0JM}{RNA levels for each gene}
#'     \item{TCGA-AO-A126}{RNA levels for each gene}
#'     \item{TCGA-AO-A12B}{RNA levels for each gene}
#'     \item{TCGA-AO-A12D}{RNA levels for each gene}
#'     \item{TCGA-AO-A12E}{RNA levels for each gene}
#'     \item{TCGA-AO-A12F}{RNA levels for each gene}
#'     \item{TCGA-AR-A0TR}{RNA levels for each gene}
#'     \item{TCGA-AR-A0TT}{RNA levels for each gene}
#'     \item{TCGA-AR-A0TV}{RNA levels for each gene}
#'     \item{TCGA-AR-A0TX}{RNA levels for each gene}
#'     \item{TCGA-AR-A0U4}{RNA levels for each gene}
#'     \item{TCGA-AR-A1AP}{RNA levels for each gene}
#'     \item{TCGA-AR-A1AS}{RNA levels for each gene}
#'     \item{TCGA-AR-A1AV}{RNA levels for each gene}
#'     \item{TCGA-AR-A1AW}{RNA levels for each gene}
#'     \item{TCGA-BH-A0AV}{RNA levels for each gene}
#'     \item{TCGA-BH-A0BV}{RNA levels for each gene}
#'     \item{TCGA-BH-A0C1}{RNA levels for each gene}
#'     \item{TCGA-BH-A0C7}{RNA levels for each gene}
#'     \item{TCGA-BH-A0DD}{RNA levels for each gene}
#'     \item{TCGA-BH-A0DG}{RNA levels for each gene}
#'     \item{TCGA-BH-A0E1}{RNA levels for each gene}
#'     \item{TCGA-BH-A0E9}{RNA levels for each gene}
#'     \item{TCGA-BH-A18N}{RNA levels for each gene}
#'     \item{TCGA-BH-A18Q}{RNA levels for each gene}
#'     \item{TCGA-BH-A18U}{RNA levels for each gene}
#'     \item{TCGA-C8-A12L}{RNA levels for each gene}
#'     \item{TCGA-C8-A12T}{RNA levels for each gene}
#'     \item{TCGA-C8-A12U}{RNA levels for each gene}
#'     \item{TCGA-C8-A12V}{RNA levels for each gene}
#'     \item{TCGA-C8-A12Z}{RNA levels for each gene}
#'     \item{TCGA-C8-A130}{RNA levels for each gene}
#'     \item{TCGA-C8-A131}{RNA levels for each gene}
#'     \item{TCGA-C8-A134}{RNA levels for each gene}
#'     \item{TCGA-C8-A135}{RNA levels for each gene}
#'     \item{TCGA-C8-A138}{RNA levels for each gene}
#'     \item{TCGA-D8-A142}{RNA levels for each gene}
#'     \item{TCGA-E2-A154}{RNA levels for each gene}
#'     \item{TCGA-E2-A158}{RNA levels for each gene}
#' }
#' @source \url{https://cptac-data-portal.georgetown.edu/cptac/s/S029}
"sample_rnadata"
