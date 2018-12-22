#'Decontaminate metabarcoding data based on optimal regressions
#'
#'Takes a data frame of metabacroding reads (structured as a column of OTU IDs, followed by at least one column of reads from blanks, followed by columns of reads from samples, optionally followed by a column of taxonomic information). Ii identifies the reads that are from contamination, then removes them. It estimates the number of overlapping OTUs between the sample and the blank, and it chooses the best equation based on that. Default (recommended) is to use the decon.primary() function anytime that the estimated overlap is <40 or >400, and to use the decon.secondary() function anytime that it is >=40 or <= 400.
#'
#'@param data A data frame of metabarcoding read data consisting of at least 3 columns in this order: a column of unique OTU names/labels, at least one column of read data from a blank sample (this contains your known contaminant reads), at least one column of read data for an actual sample (each column is a sample, each row is an OTU, and each cell is the number of reads). It can optionally include a final column with taxonomy information. If multiple blanks are included (recommended), they must be in consecutive columns, starting with column 2.
#'@param numb.blanks Numeric (default = 1). Specifies the number of blanks included in the data set (if multiple blanks are included, they must be in consecutive columns, starting with column 2).
#'@param taxa Logical (T/F). Specifies whether or not the last column contains taxonomic information (default = T)
#'@param runs Numeric (default = 2). Specifies the number of times that the function should run the decontamination procedure on the data. Based on simulation results, using two runs is best on average, but using one run is better if there is very little contamination, and using more than two runs is better if there is substantial contamination (see User's Guide section 1.4.3)
#'@param low.threshold Numeric (default = 40). Selects the lower point for switching between the two decon functions (decon.regress1() and decon.regress2()). It uses decon.regress2 anytime that the estimated overlap is <low.threshold or >up.threshold. It is usually best not to change this value.
#'@param up.threshold Numeric (default = 400). Selects the higher point for switching between the two decon functions (decon.regress1() and decon.regress2()). It uses decon.regress2 anytime that the estimated overlap is <low.threshold or >up.threshold. It is usually best not to change this value.
#'
#'@return A data frame structured the same as the input but with the contaminant reads removed (rows may have been re-ordered). Additionally, all blank columns are condensed into a single mean blank column (the mean of the blanks).
#'
#'@export decon
#'
#'
decon <- function(data,numb.blanks = 1,taxa=T,runs=2,low.threshold=40,up.threshold=400){

  if(numb.blanks >1){
    #calculates mean blank
    mean.blank <- rowMeans(data[,2:(numb.blanks+1)])
    #replaces first blank with mean blank
    data[,2] <- mean.blank
    #removes other blanks
    data <- data[,-(3:(numb.blanks+1))]
    col.names.data <- colnames(data)
    col.names.data[2] <- "Mean blank"
    colnames(data) <- col.names.data}

  #saves taxa column and column names
  column.names.save <- colnames(data)
  if(taxa == T){data.taxa <- as.data.frame(cbind(as.character(data[,1]),as.character(data[,ncol(data)])))
  colnames(data.taxa) <- c("OTU","taxa")}
  if(taxa == T){data <- as.data.frame(cbind(data[,1:(ncol(data)-1)]))}else{data <- data}

  #estimates overlap for each otu
  numb.overlap <- est.overlap(data,taxa=F)

  dec <- function(x){
    for(f in 1:(ncol(x)-2)){
      x.f <- data.frame(x[,1:2],x[,(f+2)])
      if(numb.overlap[f] < low.threshold |numb.overlap[f] >up.threshold){
        res.f <- decon.regress2(x.f,1,F)}else{
          res.f <- decon.regress1(x.f,1,F)}
      if(f == 1){results.f <- res.f[1:2]}
      res.f <- res.f[,-2]
      colnames(res.f) <- c(colnames(res.f[1]),f)
      results.f <- merge(results.f,res.f,colnames(res.f[1]),all.x =T, all.y = T)}
    results.f}


  results <- data
  for(i in 1:runs){
   results <- dec(results)
    numb.overlap <- est.overlap(results,taxa=F)}

  colnames(results) <- c("OTU",c(1:(ncol(results)-1)))

  if(taxa ==T){
    results <- merge(results,data.taxa,"OTU",all.x=T,all.y=T)
    colnames(results) <- column.names.save}else{results <- results
    colnames(results) <- column.names.save}
  results}
