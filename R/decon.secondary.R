#decontaminates sample based on regression 2
#Takes a data frame of only three columns (ID, blank, sample)

decon.regress2<- function(contamination){

  colnames(contamination)[1:2] <- c("OTU","blank")
 
  #calcualtes the proportions for each sample and the blank (for all data)
  prop <- cbind.data.frame(contamination[,1],apply(contamination[,2:ncol(contamination)],2,prop.trans))
  colnames(prop) <- c("OTU","blank","sample")

  #this subsets the proptions to just the otus that amplified in the bale
  cont.prop <- subset(prop, blank   > 0)
  cont.prop <- cont.prop[,2:ncol(cont.prop)]


  #calcualtes the percent differences between each sample and the blank
  perc.dif <- (cont.prop[,1]-cont.prop)/cont.prop[,1]

  #removes the OTU column and limits the sample to just the ones that amplified in the blank
  cont <- subset(contamination, blank   > 0)

  #This removes the lables
  cont.unlabeled <- cont[,-1]

  #this preserves the lables
  labels <- cont[,1]

  prop.trans3 <- function(x){x/(sum(x)+(0.3*sum(x)))}
  
  if(sum(cont.unlabeled[,2])>0){

  #this calcuates proportions with the denominator as sum+sum*0.3. the results will be used for figuring out the number of overlaping otuzeros
  cont.prop2 <- cbind.data.frame(prop.trans(cont.unlabeled[,1]),prop.trans3(cont.unlabeled[,2]))

  perc.dif2 <- (cont.prop2[,1]-cont.prop2)/cont.prop2[,1]

    #this is for calculating the number of otus with a perc dif >0.1 using the standard algorithms (it is intended for finding otu zero)
    sample.labeled.i <- cbind.data.frame(labels,perc.dif[,2])
    colnames(sample.labeled.i) <- c("OTU","perc_dif")
    sample.labeled.i <- sample.labeled.i[order(sample.labeled.i$perc_dif,decreasing=T),]
    sample.labeled.i <- sample.labeled.i[sample.labeled.i$perc_dif != "NaN",]
    sample.labeled.i <- sample.labeled.i[sample.labeled.i$perc_dif != 1,]
    sample.labeled.i <- sample.labeled.i[sample.labeled.i$perc_dif != -Inf,]
    r <- sample.labeled.i[,2]
    r <- rank(r)
    perc.diff.i <- cbind.data.frame(r,sample.labeled.i)
    perc.diff.i <- perc.diff.i[order(perc.diff.i$r,decreasing=T),]
    perc.diff.i <- perc.diff.i[,2:3]
    otuzero.i <- sample.labeled.i[,2]
    otuzero.i <- length(otuzero.i[otuzero.i > 0.1])


    #this uses the proportions from cont.prop2 and is intended to find the number of overlaping otus and use that to select the correct regression for the data
    sample.labeled.c <- cbind.data.frame(labels,perc.dif2[,2])
    colnames(sample.labeled.c) <- c("OTU","perc_dif")
    sample.labeled.c <- sample.labeled.c[order(sample.labeled.c$perc_dif,decreasing=T),]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != "NaN",]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != 1,]
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif != -Inf,]
    sample.labeled.c <- cbind.data.frame(sample.labeled.c[,1],sample.labeled.c[,2])
    colnames(sample.labeled.c) <- c("OTU","perc_dif")
    sample.labeled.c <- sample.labeled.c[sample.labeled.c$perc_dif > 0,]
    cont.c <- merge(cont,sample.labeled.c,"OTU",all.x=T,all.y=T)
    cont.c[is.na(cont.c)] <- 0
    cont.c <- cont.c[cont.c$perc_dif > 0,]
    if(nrow(cont.c) > 1){
      cont.unlabeled.c <- cont.c[,2:(ncol(cont.c)-1)]
      cont.prop.c <- data.frame(prop.trans(cont.unlabeled.c[,1]),prop.trans3(cont.unlabeled.c[,2]))
      perc.dif.c <- (cont.prop.c[,1]-cont.prop.c)/cont.prop.c[,1]
      sample.labeled.d <- cbind.data.frame(cont.c[,1],perc.dif.c[,2])
      colnames(sample.labeled.d) <- c("OTU","perc_dif")
      sample.labeled.d <- sample.labeled.d[order(sample.labeled.d$perc_dif,decreasing=T),]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != "NaN",]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != 1,]
      sample.labeled.d <- sample.labeled.d[sample.labeled.d$perc_dif != -Inf,]
      sample.labeled.d <- cbind.data.frame(sample.labeled.d[,1],sample.labeled.d[,2])
      colnames(sample.labeled.d) <- c("OTU","perc_dif")
      numb.overlap <- nrow(sample.labeled.d[sample.labeled.d$perc_dif > 0,])}else{numb.overlap <- 0}
    numb.overlap <- nrow(cont)-numb.overlap
    numb.overlap <- (numb.overlap-17.99)/.8246

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
      if(otuzero3.i <= 0){otuzero3.i <- 1}
      if(otuzero3.i > nrow(perc.diff.i)){otuzero3.i <- nrow(perc.diff.i)}
      if(otuzero3.i <= 0){otuzero3.i <- 1}
      otuzero3.i <- perc.diff.i[otuzero3.i,1]
      otuzero3.blank.i <- cont[cont$OTU == otuzero3.i,2]
      otuzero3.blank.ratios.i <- cont[,2]/otuzero3.blank.i
      otuzero3.correction.i <- cont[cont$OTU == otuzero3.i,3]
      otuzero3.correction.i <- otuzero3.blank.ratios.i*otuzero3.correction.i}else{otuzero3.correction.i <- cont.unlabeled[,2]}
  

  corrected <- cont[,3]-otuzero3.correction.i
  corrected[corrected < 0] <- 0
  corrected <- round(corrected)

  corrected <- cbind.data.frame(labels,cont[,2],corrected)

  colnames(corrected) <- colnames(cont)
  corrected <- rbind.data.frame(corrected,subset(contamination,blank==0))

  colnames(corrected) <- c("OTU",c(1:(ncol(corrected)-1)))
  
  corrected}