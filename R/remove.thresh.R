#' Set OTUs to zero based on thresholds
#'
#' Takes a data frame of metabarcoding reads that has gone through remove.cont(), and it finds unabundant OTUs that should be removed and sets them to zero. It only does this for OTUs that amplified in the blank sample. It applies user-specified thresholds to do this, and it does this separately for user-specified groups of individuals. Two threshold arguments are available: thresh and prop.thresh. When both options are utilized, it filters by thresh then by prop.thresh. Both arguments are run independently. Thus, an OTU will be set to zero if either condition is met, rather than if both conditions are met.
#'
#'@param data A data frame of metabarcoding read data consisting of at least 3 columns in this order: a column of unique OTU names/labels, at least one column of read data from a blank sample (this contains your known contaminant reads), at least one column of read data for an actual sample (each column is a sample, each row is an OTU, and each cell is the number of reads). It can optionally include a final column with taxonomy information. If multiple blanks are included (recommended), they must be in consecutive columns, starting with column 2. Individuals must be ordered by group (e.g., species, populations, etc.)
#'@param taxa Logical (T/F). Specifies whether or not the last column contains taxonomic information (default = T)
#'@param numb.ind A vector of numbers listing the number of individuals in each user-specified group (e.g., different populations or different species could be treated as different groups). Data must be sorted by these groups beforehand.
#'@param thresh Numeric (default = 0.7). A number written as a proportion. This is the threshold at which if that proportion of 0s are present for an OTU within a group, all samples will be set to 0 for that OTU for that group (e.g., if thresh = 0.7, then if, for a particular OTU, 70 percent of samples are 0 within a group, all samples become 0 for that OTU). The threshold always rounds down to calculate the maximum number of zeros that can be present (e.g., if thresh = 0.7 and there are 11 samples, then any OTU with 7 or more 0s will become 0 for all samples in that group). It will not do anything to groups with four or fewer samples. Set to 1 if you do not want to apply this threshold
#'@param prop.thresh Numeric (default = 0.00005). A number written as a proportion. This is the threshold at which if the number of reads for a particular OTU are below this proportion, the OTU will be set to zero for all individuals in that group (e.g., if a particular OTU makes up 0.001 percent of all of the reads for a group, then at prop.thresh = 0.00005, that OTU would be set to 0 for all individuals in the group [0.00005 = 0.005 percent]). The proportions are based on all reads for all individuals in a group (including OTUs that were not in the blank). It is necessary to relax this threshold (e.g., 0.0005) for very small data sets (see User Guide section 1.4.4) Set to 0 if you do not want to use this threshold.
#'
#'@return A data frame structured the same as the input but with the unabundant contaminant reads set to zero (rows may have been re-ordered)
#'
#'@export remove.thresh
#'
#'
remove.thresh <- function(data,numb.ind,taxa=T,thresh = 0.7,prop.thresh = 0.00005){

  if(length(unique(data[,1]))<length(data[,1])){stop("OTU ID column (column 1) contains duplicate values. All OTUs must have a unique ID")}
  if(sum(2,numb.ind,if(taxa==T){1}) != ncol(data)){stop("Total number of input columns does not match the number of columns in 'data.' Check that an OTU ID column is present, there is a single column for the blank, the number of individuals (numb.ind) are entered correctly, and the taxa column is correctly specified as present or absent")} 

  column.names <- colnames(data)
  if(taxa == T){data.taxa <- cbind.data.frame(data[,1],data[,ncol(data)])
  colnames(data.taxa) <- c("OTU","taxa")}
  if(taxa == T){data <- cbind.data.frame(data[,1:(ncol(data)-1)])}

  header <- rep("samp",ncol(data)-2)
  header <- gsub(" ","",paste(header,1:(ncol(data)-2)))
  header <- c("OTU","blank",header)
  colnames(data) <- header

  #this calculates the total number of reads for each group (before anything is removed)
  group.sum <- NULL
  for(c in 1:length(numb.ind)){
    group.c <- data[,(sum(numb.ind[1:c])-numb.ind[c]+3):(sum(numb.ind[1:c])+2)]
    group.c <- sum(group.c)
    group.sum <- c(group.sum,group.c)}


  cont <- data[data$blank > 0,]


  count.0 <- function(x){length(x[x==0])}


  corrected <- as.data.frame(cont[,3:ncol(cont)])
  corrected.a <- cont[,1:2]
  for(a in 1:length(numb.ind)){
    pop.a <- cbind.data.frame(corrected[,(sum(numb.ind[1:a])-numb.ind[a]+1):(sum(numb.ind[1:a]))])

    #this condition prevents it from making corrections based on the number of zeros (thresh) if fewer than four samples are present, but it will still correct based on prop.thresh
    if(ncol(pop.a) > 4){
      pop.a <- cbind.data.frame(pop.a,apply(pop.a,1,count.0))

      pop.b <- NULL
      for(b in 1:nrow(pop.a)){
        row.b <- pop.a[b,]
        if(pop.a[b,ncol(pop.a)] >= floor(thresh*(ncol(pop.a)-1)) | ((sum(pop.a[b,1:(ncol(pop.a)-1)]))/group.sum[a]) < prop.thresh){row.b[row.b > 0] <- 0}else{row.b <- row.b}
        pop.b <- rbind.data.frame(pop.b,row.b)}
     pop.a <- pop.b[,-ncol(pop.b)]}else{pop.a <- pop.a}

    pop.b <- NULL
    for(b in 1:nrow(pop.a)){
      row.b <- pop.a[b,]
      if(((sum(pop.a[b,1:(ncol(pop.a))]))/group.sum[a]) < prop.thresh){row.b[row.b > 0] <- 0}else{row.b <- row.b}
      pop.b <- rbind.data.frame(pop.b,row.b)}

    corrected.a <- cbind.data.frame(corrected.a,pop.b)}

  colnames(corrected.a) <- colnames(cont)
  corrected.a <- rbind.data.frame(corrected.a,subset(data,blank==0))

  colnames(corrected.a) <- c("OTU",c(1:(ncol(corrected.a)-1)))

  if(taxa ==T){
    corrected.a <- merge(corrected.a,data.taxa,"OTU",all.x=T,all.y=T)
    colnames(corrected.a) <- column.names}else{corrected.a <- corrected.a
    colnames(corrected.a) <- column.names}
  corrected.a}