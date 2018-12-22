#' Set OTUs to zero based on threshold
#'
#' Takes a data frame of metabarcoding reads that has gone through one of the decon() functions, and it finds unabundant OTUs that should be removed and sets them to zero. It only does this for OTUs that amplified in the blank sample. It applies user-specified thresholds to do this, and it does this separately for user-specified groups of individuals. Two threshold arguments are available: thresh and prop.thresh. When both options are utilized, it filters by “thresh” then by “prop.thresh.” Both arguments are run independently. Thus, an OTU will be set to zero if either condition is met, rather than if both conditions are met.
#'
#'@param data A data frame of metabarcoding read data consisting of at least 3 columns in this order: a column of unique OTU names/labels, a column of read data from a blank sample (this contains your known contaminant reads; do not include multiple blank columns; include only the mean blank column), at least one column of read data for an actual sample (each column is a sample, each row is an OTU, and each cell is the number of reads). It can optionally include a final column with taxonomy information. If multiple groups are present, order all your sample columns so that groups are together, and include the corresponding number of individuals per group in numb.ind.
#'@param taxa Logical (T/F). Specifies whether or not the last column contains taxonomic information (default = T)
#'@param numb.ind A vector of numbers listing the number of individuals in each user-specified group (e.g., different populations could be treated as different groups). Data must be sorted by these groups beforehand. To apply the thresholds across all individuals as a single group, make a vector with a single number that is the total number of individuals. Alternatively, to treat each sample independently, you can set each individual as a group. Do this by making a vector with the number 1 repeated as many times as you have samples (e.g., if you have 20 samples, you can make this as follows: my.vector <- rep(1,20)). If each individual is treated separately, then only the prop.thresh argument will be used for filtering.
#'@param thresh Numeric (default = 0.7). A number written as a proportion. This is the threshold at which if that proportion of 0s are present for an OTU within a group, all samples will be set to 0 for that group (e.g., if thresh = 0.7, then if, for a particular OTU, 70 percent of samples are 0 within a group, all become 0 for that OTU). The threshold always rounds down (i.e. if thresh = 0.7 and there are 11 samples, then any OTU with 7 or more 0s will become 0 for all samples in that group). It will not do anything to groups with four or fewer samples. Set to 1 if you don't want to apply this threshold.
#'@param prop.thresh Numeric (default = 0.00005). A number written as a proportion. This is the threshold at which if the number of reads for a particular OTU are below this proportion (based on all reads for a group), the OTU will be set to zero for all individuals in that group (e.g., if a particular OTU only makes up 0.001 percent of all of the reads for a group, then at prop.thresh = 0.00005, that OTU would be set to 0 for all individuals in the group [0.00005 = 0.005 percent]). Set to 0 if you don't want to use this threshold.
#'
#'@return A data frame structured the same as the input but with the unabundant contaminant reads set to zero (rows may have been re-ordered).
#'
#'@export remove.thresh
#'
remove.thresh <- function(data,taxa=T,numb.ind,thresh = 0.7,prop.thresh = 0.00005){

  column.names <- colnames(data)
  if(taxa == T){data.taxa <- as.data.frame(cbind(as.character(data[,1]),as.character(data[,ncol(data)])))
  colnames(data.taxa) <- c("OTU","taxa")}
  if(taxa == T){data <- as.data.frame(cbind(data[,1:(ncol(data)-1)]))}else{data <- data}

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
    pop.a <- as.data.frame(cbind(corrected[,(sum(numb.ind[1:a])-numb.ind[a]+1):(sum(numb.ind[1:a]))]))

    #this condition prevents it from making corrections based on the number of zeros (thresh) if fewer than four samples are present, but it will still correct based on prop.thresh
    if(ncol(pop.a) > 4){
      pop.a <- as.data.frame(cbind(pop.a,apply(pop.a,1,count.0)))

      pop.b <- NULL
      for(b in 1:nrow(pop.a)){
        row.b <- pop.a[b,]
        if(pop.a[b,ncol(pop.a)] >= floor(thresh*(ncol(pop.a)-1)) | ((sum(pop.a[b,1:(ncol(pop.a)-1)]))/group.sum[a]) < prop.thresh){row.b[row.b > 0] <- 0}else{row.b <- row.b}
        pop.b <- as.data.frame(rbind(pop.b,row.b))}
     pop.a <- pop.b[,-ncol(pop.b)]}else{pop.a <- pop.a}

    pop.b <- NULL
    for(b in 1:nrow(pop.a)){
      row.b <- pop.a[b,]
      if(((sum(pop.a[b,1:(ncol(pop.a))]))/group.sum[a]) < prop.thresh){row.b[row.b > 0] <- 0}else{row.b <- row.b}
      pop.b <- as.data.frame(rbind(pop.b,row.b))}

    corrected.a <- as.data.frame(cbind(corrected.a,pop.b))}

  colnames(corrected.a) <- colnames(cont)
  corrected.a <- rbind(corrected.a,subset(data,blank==0))

  colnames(corrected.a) <- c("OTU",c(1:(ncol(corrected.a)-1)))

  if(taxa ==T){
    corrected.a <- merge(corrected.a,data.taxa,"OTU",all.x=T,all.y=T)
    colnames(corrected.a) <- column.names}else{corrected.a <- corrected.a
    colnames(corrected.a) <- column.names}
  corrected.a}



