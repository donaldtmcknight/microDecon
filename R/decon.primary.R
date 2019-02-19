#'Decontaminate metabarcoding data based on a single regression
#'
#'Takes a data frame of metabarcoding reads, identifies the reads that are from contamination, then removes them. It does this based on a single regression line that is used for each sample. It is generally recommended to use decon() rather than decon.regress1()
#'
#'@param data A data frame of metabarcoding read data consisting of at least 3 columns in this order: a column of unique OTU names/labels, at least one column of read data from a blank sample (this contains your known contaminant reads), at least one column of read data for an actual sample (each column is a sample, each row is an OTU, and each cell is the number of reads). It can optionally include a final column with taxonomy information. If multiple blanks are included (recommended), they must be in consecutive columns, starting with column 2.
#'@param numb.blanks Numeric (default = 1). Specifies the number of blanks included in the data set (if multiple blanks are included, they must be in consecutive columns, starting with column 2).
#'@param taxa Logical (T/F). Specifies whether or not the last column contains taxonomic information (default = T)
#'
#'@return A data frame structured the same as the input but with the contaminant reads removed (rows may have been re-ordered).Additionally, all blank columns are condensed into a single mean blank column (the mean of the blanks).
#'
#'@export decon.regress1
#'
#'
decon.regress1<- function(data,numb.blanks = 1,taxa=T){

  if(numb.blanks >1){
    mean.reads.blank <- mean(colSums(data[,2:(numb.blanks+1)]))
    #calculates mean blank
    mean.blank <- rowMeans(apply(data[,2:(numb.blanks+1)],2,prop.trans))*mean.reads.blank
    #replaces first blank with mean blank
    data[,2] <- mean.blank
    #removes other blanks
    data <- data[,-(3:(numb.blanks+1))]
    col.names.data <- colnames(data)
    col.names.data[2] <- "Mean blank"
    colnames(data) <- col.names.data}

  column.names <- colnames(data)
  if(taxa == T){data.taxa <- as.data.frame(cbind(as.character(data[,1]),as.character(data[,ncol(data)])))
  colnames(data.taxa) <- c("OTU","taxa")}
  if(taxa == T){data <- as.data.frame(cbind(data[,1:(ncol(data)-1)]))}else{data <- data}

  header <- rep("samp",ncol(data)-2)
  header <- gsub(" ","",paste(header,1:(ncol(data)-2)))
  header <- c("OTU","blank",header)
  colnames(data) <- header

  #removes the OTU column and limits the sample to just the ones that amplified in the blank
  cont <- subset(data, blank   > 0)

  #this preserves the lables
  labels <- cont[,1]

  #This removes the lables
  cont.unlabeled <- cont[,2:ncol(data)]

  prop.trans2 <- function(x){x/(sum(x)+(0.1*sum(x)))}

  #calcualtes the proportions for each sample and the blank (still for subset data)
  cont.prop <- cbind(cont.unlabeled[,1]/(sum(cont.unlabeled[,1])),apply(cbind(cont.unlabeled[,2:ncol(cont.unlabeled)]),2,prop.trans2))

  #calcualtes the percent differences between each sample and the blank
  perc.dif <- (cont.prop[,1]-cont.prop)/cont.prop[,1]

  correction.i <- NULL
  for(i in 1:(ncol(perc.dif)-1)){
    sample.labeled.i <- as.data.frame(cbind(as.character(labels),perc.dif[,i+1]))
    colnames(sample.labeled.i) <- c("OTU","perc_dif")
    sample.labeled.i <- sample.labeled.i[order(sample.labeled.i$perc_dif,decreasing=T),]
    sample.labeled.i <- sample.labeled.i[sample.labeled.i$perc_dif != "NaN",]
    sample.labeled.i <- sample.labeled.i[sample.labeled.i$perc_dif != 1,]
    sample.labeled.i <- sample.labeled.i[sample.labeled.i$perc_dif != -Inf,]
    r <- as.numeric(as.character(sample.labeled.i[,2]))
    r <- rank(r)
    perc.diff.i <- as.data.frame(cbind(r,sample.labeled.i))
    perc.diff.i <- perc.diff.i[order(perc.diff.i$r,decreasing=T),]
    perc.diff.i <- perc.diff.i[,2:3]
    otuzero.i <- as.numeric(as.character(sample.labeled.i[,2]))
    otuzero.i <- length(otuzero.i[otuzero.i > 0.1])

    if(otuzero.i > 0){

      otuzero3.i <- (0.7754*otuzero.i)-4.2185
      otuzero3.i <- round(otuzero3.i)
      if(otuzero3.i <= 0){otuzero3.i <- 1}else{otuzero3.i<-otuzero3.i}
      if(otuzero3.i > nrow(perc.diff.i)){otuzero3.i <- nrow(perc.diff.i)}else{otuzero3.i <- otuzero3.i}
      otuzero3.i <- as.character(perc.diff.i[otuzero3.i,1])
      otuzero3.blank.i <- cont[cont$OTU == otuzero3.i,2]
      otuzero3.blank.ratios.i <- cont[,2]/otuzero3.blank.i
      otuzero3.correction.i <- cont[cont$OTU == otuzero3.i,i+2]
      otuzero3.correction.i <- otuzero3.blank.ratios.i*otuzero3.correction.i

      mean.i <- otuzero3.correction.i}else{mean.i <- c(rep(0,nrow(perc.dif)))}

    correction.i <- as.data.frame(cbind(correction.i,mean.i))}

  corrected <- cont[,3:ncol(cont)]-correction.i
  corrected[corrected < 0] <- 0
  corrected <- round(corrected)


  corrected.labels <- as.data.frame(cbind(as.character(labels)))
  corrected <- cbind(corrected.labels,cont[,2],corrected)

  colnames(corrected) <- colnames(cont)
  corrected <- rbind(corrected,subset(data,blank==0))

  colnames(corrected) <- c("OTU",c(1:(ncol(corrected)-1)))

  if(taxa ==T){
    corrected <- merge(corrected,data.taxa,"OTU",all.x=T,all.y=T)
    colnames(corrected) <- column.names}else{corrected <- corrected
    colnames(corrected) <- column.names}
  corrected}



