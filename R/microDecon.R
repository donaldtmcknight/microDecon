#'microDecon: A package for removing contaminant sequencing reads from metabarcoding studies
#'
#'The microDecon package provides four functions for removing contamination.
#'
#'decon() = The primary function of the package. It automatically selects between the algorithms in decon.regress1() and decon.regress2(). It is recommended to use this function on default settings in conjunction with the remove.thresh() function.
#'
#'decon.regress1() = A function that uses a single, fixed algorithm for removing contamination. It works well in most cases, but it is recommended to use decon() instead.
#'
#'decon.regress2() = A function that uses a single, fixed algorithm for removing contamination. It only works well when there are either very few or very many overlapping OTUs. It is recommended to use decon() instead.
#'
#'remove.thresh() = A function to remove residual contamination following the use of one of the decon functions.
#'
#'@author Donald T. McKnight <donald.mcknight@my.jcu.edu.au>
#'
#'
#'@docType package
#'@name microDecon
#'@keywords internal

NULL
