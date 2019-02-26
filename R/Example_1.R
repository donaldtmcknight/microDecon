#' An example contaminated OTU table 
#'
#' An OTU table containing data for three blanks, three samples (individuals from a single population), and taxa, formatted for microDecon. OTU2 is present as both contamination and actual reads, whereas OTUs 3, 4, and 6 are entirely contamination
#'
#' @format A data frame with 6 rows and 8 variables:
#' \describe{
#'   \item{OTU_ID}{column of unique OTU IDs}
#'   \item{Blank1}{reads from Blank 1}
#'   \item{Blank2}{reads from Blank 2}
#'   \item{Blank3}{reads from Blank 3}
#'   \item{Ind1}{reads from individual 1}
#'   \item{Ind2}{reads from individual 2}
#'   \item{Ind3}{reads from individual 3}
#'   \item{taxa}{taxonomic information}
#' }
#' @source in silico data
"Example_1"