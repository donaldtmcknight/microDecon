#'Calculates summary statistics for results of decontamination
#'
#'Takes your original data and decontaminated data and returns summary statistics such as the reads removed for each OTU and sample ($reads.removed), mean reads removed for each group and all samples ($difference.mean), sum of reads removed for each group and all samples ($difference.sum), IDs of OTUs that were totally removed from at least one group ($OTUs.removed), and the decontaminated data frame after removing any OTUs that were entirely contamination ($decon.table).
#'
#'@param data The original data frame that you input into remove.cont(). Individuals must have been ordered by groups (populations, species, etc.) as in remove.thresh().
#'@param output The data frame that was returned by remove.cont() or remove.thresh().
#'@param numb.blanks Numeric (default = 1). Specifies the number of blanks included in the “data” argument (if multiple blanks are included, they must be in consecutive columns, starting with column 2). This only applies to the number in the data argument. The number in the output argument will always be 1 because remove.cont() returns a single mean blank.
#'@param numb.ind A vector of numbers listing the number of individuals in each user-specified group (e.g., different populations could be treated as different groups). Data must have been sorted by these groups before running remove.cont().
#'@param taxa Logical (T/F). Specifies whether or not the last column contains taxonomic information (default = T).
#'
#'@return A list of five data frames that can be accessed with $. These are useful for both seeing and recording the changes microDecon made, as well as checking that the changes make sense based on the biological understanding of the system under study.
#'
#'  NA values indicate that an OTU had zero reads for a given group or sample prior to decontamination.
#'  
#'  $decon.table = A data frame of decontaminated OTU data. It is structured the same as the original data frame (data). However, if several blanks were input, the output will include only a single Mean.blank column that is the mean of the proportions of those blanks multiplied by the mean number of reads in the blanks. Additionally, the order of the rows may be different, and any OTUs for which all reads were removed will have been deleted (their information will still be shown in the other outputs).
#'  
#'  $reads.removed = An OTU table showing the number of reads that were removed from each OTU that amplified in the blank (per individual).
#'  
#'  $difference.sum = The total number of reads that were removed from each OTU that amplified in the blank (per group as well as for the entire data set; groups are in the same order as specified by the numb.ind argument).
#'  
#'  $difference.mean = The average number of reads that were removed from each OTU that amplified in the blank (per group as well as for the entire data set; groups are in the same order as specified by the numb.ind argument).
#'  
#'  $OTUs.removed = A data frame showing the identities of OTUs that were completely removed from either particular groups or the entire data set.
#'  
#'@export decon.diff
#'
#'
decon.diff <- function(data,output,numb.blanks=1,numb.ind,taxa=T){
 
  if(sum(1,numb.blanks,numb.ind,if(taxa==T){1}) != ncol(data)){stop("Total number of input columns does not match the number of columns in 'data.' Check that an OTU ID column is present, the number of blanks (numb.blanks) and number of individuals (numb.ind) are entered correctly, and the taxa column is correctly specified as present or absent")} 

  if(numb.blanks >1){
    #calcualtes mean reads per blank
    mean.reads.blank <- mean(colSums(data[,2:(numb.blanks+1)]))
    #calculates mean blank
    mean.blank <- rowMeans(apply(data[,2:(numb.blanks+1)],2,prop.trans))*mean.reads.blank
    #replaces first blank with mean blank
    data[,2] <- mean.blank
    #removes other blanks
    data <- data[,-(3:(numb.blanks+1))]
    col.names.data <- colnames(data)
    col.names.data[2] <- "Mean.blank"
    colnames(data) <- col.names.data}
  
  cnames <- colnames(data)[3:ncol(data)]

  data <- data[data[,2] > 0,] #subsets to data in blank
  output.sub <- output[output[,2] > 0,] #subsets to data in blank

  data <- data[order(data[,1]),] #orders data
  output.sub <- output.sub[order(output.sub[,1]),] #orders data
  
  if(taxa==T){save.taxa <- cbind.data.frame(data[,1],data[,ncol(data)]) #saves taxa
    colnames(save.taxa) <- c(colnames(data)[1],colnames(data)[ncol(data)])
    data <- data[, -ncol(data)] #removes taxa
    output.sub <- output.sub[,-ncol(output.sub)]} #removes taxa
  
  data.sub <- data[,3:ncol(data)]
  data.sub[data.sub == 0] <- NA
  data <- cbind.data.frame(data[,1:2],data.sub)
  
  difference <- cbind.data.frame(data[,1:2],(data[,3:ncol(data)]-output.sub[,3:ncol(output.sub)])) #calcualtes number of reads removed

  if(ncol(difference) > 3){
    difference.sum <- cbind.data.frame(difference[,1:2],rowSums(difference[,3:ncol(difference)],na.rm=T))
    difference.mean <- cbind.data.frame(difference[,1:2],rowMeans(difference[,3:ncol(difference)],na.rm=T))}else{
    difference.sum <- difference
      difference.mean <- difference}  
 
  removed <- cbind.data.frame(data[,3:ncol(data)]-difference[,3:ncol(difference)]) # subtracts reads removed from data and sums (0=all reads removed)
 
  l.0 <- function(x){length((x[x==0]))} #function for counting the number of zeros presant
  l.NA <- function(x){length(x[(is.na(x)==T)])} #function for counting the number of NA presant (these were 0 prior to decon)
  true.0 <- function(x){if(x[1] > 0 && x[2]==1){1}else{0}} #checks for whether all reads were removed for a given OTU for a data frame where the first column is number of individuals that were totally removed, and the second is the number of individuals and NA values
  
  total.removed <- cbind.data.frame(apply(removed,1,l.0)-apply(removed,1,l.NA),apply(removed,1,l.0)/ncol(removed)) #counts number of zeros from decon and total number of zeros (NA) (from original data and decon) divided by number of individals
  total.removed <- apply(total.removed,1,true.0)
  
  #for each group returns a 0 if all OTUs were removed, combines with total.remvoed
  for(g in 1:length(numb.ind)){
    if(g==1){sum.g <- 1} #sum.g = sum of previous numbers of individuals (set to 1 for round 1)
    group.g <- cbind.data.frame(removed[,sum.g:(sum.g+numb.ind[g]-1)])
    group.g <- cbind.data.frame(apply(group.g,1,l.0)-apply(group.g,1,l.NA),apply(group.g,1,l.0)/ncol(group.g)) #counts number of zeros from decon and total number of zeros (NA) (from original data and decon) divided by number of individals
    group.g <- apply(group.g,1,true.0)
    total.removed <- cbind.data.frame(total.removed,group.g)
    
    if(numb.ind[g] > 1){diff.sum.g <- rowSums(difference[,(sum.g+2):(sum.g+numb.ind[g]+1)],na.rm=T)}else{diff.sum.g <- difference[,(sum.g+2)]}
    difference.sum <- cbind.data.frame(difference.sum,diff.sum.g)
    
    diff.mean.g <- rowMeans(cbind.data.frame(difference[,(sum.g+2):(sum.g+numb.ind[g]+1)]),na.rm=T)
    difference.mean <- cbind.data.frame(difference.mean,diff.mean.g)
    
    sum.g <- sum.g+numb.ind[g]}
  
  difference.mean[is.na(difference.mean)==T] <- NA #relaces nan values from otus that were not present to beign with
  
  total.removed <- cbind.data.frame(data[,1:2],total.removed,rowSums(total.removed)) #adds column with number of groups from which it was totally removed
  colnames(total.removed)[ncol(total.removed)] <- "rem"
total.removed <- total.removed[total.removed$rem > 0,] #subsets to OTUs where something had been totally removed at some point
  if(nrow(total.removed) > 0){
  
    total.removed <-  total.removed[,-ncol( total.removed )] #removes column with info on number totally removed
    
    #orders total.removed to be the same as difference.mean, then uses the NA locations in difference.mean to assign NA to total.removed
    total.removed <- total.removed[order(total.removed[,1]),]
    diff.na <- subset(  difference.mean,difference.mean[,1]%in% total.removed[,1])
    diff.na <- diff.na[order(diff.na[,1]),]
    total.removed[is.na(diff.na)==T] <- NA
    
    total.removed.otu <- total.removed[,1:2] #saves otu and blank column for subset data
    total.removed <-total.removed[,3:ncol(total.removed)]
  
    total.removed[total.removed == 0] <- "-"
    total.removed[total.removed == 1] <- "Totally.removed"
    total.removed <- cbind.data.frame(total.removed.otu,total.removed)
    
    colnames(total.removed)[3:ncol(total.removed)] <- c("All.groups",gsub(" ","",paste("Group",c(1:length(numb.ind)))))
    
    if(taxa == T){ #adds taxa info
      total.removed <- merge(total.removed,save.taxa,colnames(save.taxa)[1],all.x=F,all.y=F)}
    
     result.removed <- setdiff(output[,1],total.removed[total.removed$All.groups == "Totally.removed",1]) #returns vector of OTUs that were not totally removed
     result.removed <- subset(output,output[,1]%in%result.removed)}else{ #subsets data to those OTUs
  
      result.removed <- output
      total.removed <- "no OTUs were totally removed from any groups"}
  
  if(taxa == T){
    colnames(difference.mean)[3:ncol(difference.mean)] <- c("All.groups",gsub(" ","",paste("Group",c(1:length(numb.ind)))))
    colnames(difference.sum)[3:ncol(difference.sum)] <- c("All.groups",gsub(" ","",paste("Group",c(1:length(numb.ind)))))
    
    difference.mean <- merge(difference.mean,save.taxa,by=colnames(save.taxa)[1])
    difference.sum <- merge(difference.sum,save.taxa,by=colnames(save.taxa)[1])
    difference <- merge(difference,save.taxa,by=colnames(save.taxa)[1])
  
    }else{
      colnames(difference.mean)[3:ncol(difference.mean)] <- c("All.groups",gsub(" ","",paste("Group",c(1:length(numb.ind)))))
      colnames(difference.sum)[3:ncol(difference.sum)] <- c("All.groups",gsub(" ","",paste("Group",c(1:length(numb.ind)))))}     
  
  colnames(difference)[3:ncol(difference)] <- cnames
  
  results <- list("reads.removed"=difference,"mean.per.group"=difference.mean,"sum.per.group"=difference.sum,"OTUs.removed"=total.removed,"decon.table" = result.removed)}
  