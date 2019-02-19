#' Decontaminate metabarcoding data based on several regressions
#'
#' Takes a data frame of metabarcoding reads, identifies the reads that are from contamination, then removes them. It does this based on several different regression lines. It chooses which line to use for a particular sample based on its estimate of the number of overlapping OTUs between the blank and the sample. It is not recommended except possibly in cases where there are very few overlapping OTUs or where there are >400 overlapping OTUs (generally use the decon() function).
#'
#'@param data A data frame of metabarcoding read data consisting of at least 3 columns in this order: a column of unique OTU names/labels, at least one column of read data from a blank sample (this contains your known contaminant reads), at least one column of read data for an actual sample (each column is a sample, each row is an OTU, and each cell is the number of reads). It can optionally include a final column with taxonomy information. If multiple blanks are included (recommended), they must be in consecutive columns, starting with column 2.
#'@param numb.blanks Numeric (default = 1). Specifies the number of blanks included in the data set (if multiple blanks are included, they must be in consecutive columns, starting with column 2).
#'@param taxa Logical (T/F). Specifies whether or not the last column contains taxonomic information (default = T)
#'
#'@return A data frame structured the same as the input but with the contaminant reads removed (rows may have been re-ordered). Additionally, all blank columns are condensed into a single mean blank column (the mean of the blanks).
#'
#'@export decon.regress2
#'
#'
decon.regress2<- function(data,numb.blanks = 1,taxa=T){

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

  prop.trans <- function(x){x/sum(x)}

  #calcualtes the proportions for each sample and the blank (for all data)
  prop <- apply(data[,2:ncol(data)],2,prop.trans)
  prop <- data.frame(data[,1],prop)
  colnames(prop) <- header

  #this subsets the proptions to just the otus that amplified in the bale
  cont.prop <- subset(prop, blank   > 0)
  cont.prop <- cont.prop[,2:ncol(cont.prop)]


  #calcualtes the percent differences between each sample and the blank
  perc.dif <- (cont.prop[,1]-cont.prop)/cont.prop[,1]

  #removes the OTU column and limits the sample to just the ones that amplified in the blank
  cont <- subset(data, blank   > 0)

  #This removes the lables
  cont.unlabeled <- cont[,2:ncol(data)]

  #this preserves the lables
  labels <- cont[,1]

  prop.trans3 <- function(x){x/(sum(x)+(0.3*sum(x)))}

  #this calcuates proportions with the denominator as sum+sum*0.3. the results will be used for figuring out the number of overlaping otuzeros
  cont.prop2 <- data.frame(apply(cbind(cont.unlabeled[,1]),2,prop.trans),apply(cbind(cont.unlabeled[,2:ncol(cont.unlabeled)]),2,prop.trans3))

  perc.dif2 <- (cont.prop2[,1]-cont.prop2)/cont.prop2[,1]

  correction.i <- NULL
  for(i in 1:(ncol(perc.dif)-1)){

    #this is for calculating the number of otus with a perc dif >0.1 using the standard algorythims (it is intended for finding otu zero)
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


    #this uses the proportions from cont.prop2 and is intended to find the number of overlaping otus and use that to select the correct regression for the data
    sample.labeled.c <- as.data.frame(cbind(as.character(labels),perc.dif2[,i+1]))
    colnames(sample.labeled.c) <- c("OTU","perc_dif")
    sample.labeled.c <- sample.labeled.c[order(sample.labeled.c$perc_dif,decreasing=T),]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != "NaN",]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != 1,]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != -Inf,]
    sample.labeled.c <- data.frame(as.character(sample.labeled.c[,1]),as.numeric(as.character(sample.labeled.c[,2])))
    colnames(sample.labeled.c) <- c("OTU","perc_dif")
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif > 0,]
    cont.c <- merge(cont,sample.labeled.c,"OTU",all.x=T,all.y=T)
    cont.c[is.na(cont.c)] <- 0
    cont.c <- cont.c[cont.c$perc_dif > 0,]
    if(nrow(cont.c) > 1){
      cont.unlabeled.c <- cont.c[,2:(ncol(cont.c)-1)]
      cont.prop.c <- data.frame(apply(cbind(cont.unlabeled.c[,1]),2,prop.trans),apply(cbind(cont.unlabeled.c[,2:ncol(cont.unlabeled.c)]),2,prop.trans3))
      perc.dif.c <- (cont.prop.c[,1]-cont.prop.c)/cont.prop.c[,1]
      sample.labeled.d <- data.frame(cont.c[,1],perc.dif.c[,i+1])
      colnames(sample.labeled.d) <- c("OTU","perc_dif")
      sample.labeled.d <- sample.labeled.d[order(sample.labeled.d$perc_dif,decreasing=T),]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != "NaN",]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != 1,]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != -Inf,]
      sample.labeled.d <- data.frame(as.character(sample.labeled.d[,1]),as.numeric(as.character(sample.labeled.d[,2])))
      colnames(sample.labeled.d) <- c("OTU","perc_dif")
      numb.overlap <- nrow(sample.labeled.d[sample.labeled.d$perc_dif > 0,])}else{numb.overlap <- 0}
    numb.overlap <- nrow(cont)-numb.overlap
    numb.overlap <- (numb.overlap-17.99)/.8246

    if(otuzero.i > 0){

      if(numb.overlap < 21){m <- 0.953
      b <- -10.715}
      if(numb.overlap > 20 & numb.overlap < 41){m <- 0.8833
      b <- -22.003}
      if(numb.overlap > 40 & numb.overlap < 61){m <- 0.896
      b <- -33.91}
      if(numb.overlap > 60 & numb.overlap < 81){m <- 1.0393
      b <- -65.656}
      if(numb.overlap > 80 & numb.overlap < 101){m <- 1.1773
      b <- -93.428}
      if(numb.overlap > 100 & numb.overlap < 151){m <- 1.1827
      b <- -126.22}
      if(numb.overlap > 150 & numb.overlap < 201){m <- 1.1608
      b <- -159.31}
      if(numb.overlap > 200){m <- 0.9764
      b <- -172.22}

      otuzero3.i <- (m*otuzero.i)+b
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



