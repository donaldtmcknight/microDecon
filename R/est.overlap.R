#takes a data frame with read data (number of reads per otu per sample) with the following columns (data should be sorted by groups prior to runing this function):
# Column 1 = OTU IDs (e.g., otu1, otu2). The exact names are irrelevant, but each row must have a unique one
# Column 2 = reads from a blank sample
# Columns 3-i = Columns of reads for each sample (each column is a sample)


#data = your data frame

est.overlap <- function(data){

  colnames(data)[1:2] <- c("OTU","blank")

  #removes the OTU column and limits the sample to just the ones that amplified in the blank
  cont <- subset(data, blank   > 0)

  #this preserves the lables
  labels <- cont[,1]

  #This removes the lables
  cont.unlabeled <- cont[,2:ncol(data)]

  prop.trans3 <- function(x){x/(sum(x)+(0.3*sum(x)))}

  #this calcuates proportions with the denominator as sum+sum*0.3. the results will be used for figuring out the number of overlaping otuzeros
  cont.prop2 <- cbind.data.frame(prop.trans(cont.unlabeled[,1]),apply(cbind(cont.unlabeled[,2:ncol(cont.unlabeled)]),2,prop.trans3))

  perc.dif2 <- (cont.prop2[,1]-cont.prop2)/cont.prop2[,1]

  
  for(i in 1:(ncol(perc.dif2)-1)){

    #this uses the proportions from cont.prop2 and is intended to find the number of overlaping otus and use that to select the correct regression for the data
    sample.labeled.c <- cbind.data.frame(labels,perc.dif2[,i+1])
    colnames(sample.labeled.c) <- c("OTU","perc_dif")
    sample.labeled.c <- sample.labeled.c[order(sample.labeled.c$perc_dif,decreasing=T),]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != "NaN",]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != 1,]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != -Inf,]
    sample.labeled.c <- sample.labeled.c[,1:2]
    colnames(sample.labeled.c) <- c("OTU","perc_dif")
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif > 0,]
    cont.c <- merge(cont,sample.labeled.c,"OTU",all.x=T,all.y=T)
    cont.c[is.na(cont.c)] <- 0
    cont.c <- cont.c[cont.c$perc_dif > 0,]
    if(nrow(cont.c) > 1){
      cont.unlabeled.c <- cont.c[,2:(ncol(cont.c)-1)]
      cont.prop.c <- cbind.data.frame(prop.trans(cont.unlabeled.c[,1]),apply(cbind(cont.unlabeled.c[,2:ncol(cont.unlabeled.c)]),2,prop.trans3))
      perc.dif.c <- (cont.prop.c[,1]-cont.prop.c)/cont.prop.c[,1]
      sample.labeled.d <- data.frame(cont.c[,1],perc.dif.c[,i+1])
      colnames(sample.labeled.d) <- c("OTU","perc_dif")
      sample.labeled.d <- sample.labeled.d[order(sample.labeled.d$perc_dif,decreasing=T),]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != "NaN",]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != 1,]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != -Inf,]
      sample.labeled.d <- sample.labeled.d[,1:2]
      colnames(sample.labeled.d) <- c("OTU","perc_dif")
      numb.overlap <- nrow(sample.labeled.d[sample.labeled.d$perc_dif > 0,])}else{numb.overlap <- nrow(cont.c)}
    numb.overlap <- nrow(cont)-numb.overlap
    numb.overlap <- (numb.overlap-17.99)/.8246
    mean.i <- numb.overlap

    correction.i <- if(i==1){mean.i}else{cbind.data.frame(correction.i,mean.i)}}

  correction.i}