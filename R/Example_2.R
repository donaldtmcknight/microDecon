#' An example contaminated OTU table 
#'
#' An OTU table containing data for one blank and 10 samples (six individuals from one population and four from another), formatted for microDecon. OTU2 is actually present on samples in group 1 (but not group 2), and OTU4 is actually present in group 2 (but not group 1). OTUs 3 and 6 are entirely contamination. 
#'
#' @format A data frame with 6 rows and 8 variables:
#' \describe{
#'   \item{OTU_ID}{column of unique OTU IDs}
#'   \item{Blank}{reads from the blank}
#'   \item{Ind1}{reads from individual 1 (population 1)}
#'   \item{Ind2}{reads from individual 2 (population 1)}
#'   \item{Ind3}{reads from individual 3 (population 1)}
#'   \item{Ind4}{reads from individual 4 (population 1)}  
#'   \item{Ind5}{reads from individual 5 (population 1)}  
#'   \item{Ind6}{reads from individual 6 (population 1)}  
#'   \item{Ind7}{reads from individual 7 (population 2)}  
#'   \item{Ind8}{reads from individual 8 (population 2)}  
#'   \item{Ind9}{reads from individual 9 (population 2}  
#'   \item{Ind10}{reads from individual 10 (population 2)}     
#' }
#' @source in silico data
"Example_2"