#decontaminates sample based on regression 1 
#Takes a data frame of only three columns (ID, blank, sample)

decon.regress1<- function(contamination){

 
 colnames(contamination)[1:2] <- c("OTU","blank")
 
  #removes the OTU column and limits the sample to just the ones that amplified in the blank
  cont <- subset(contamination, blank   > 0)

  #this preserves the lables
  labels <- cont[,1]

  #This removes the lables
  cont.unlabeled <- cont[,2:ncol(contamination)]

  if(sum(cont.unlabeled[,2])>0){
  
  prop.trans2 <- function(x){x/(sum(x)+(0.1*sum(x)))}

  #calcualtes the proportions for each sample and the blank (still for subset contamination)
  cont.prop <- cbind.data.frame(prop.trans(cont.unlabeled[,1]),prop.trans2(cont.unlabeled[,2]))

  #calcualtes the percent differences between each sample and the blank
  perc.dif <- (cont.prop[,1]-cont.prop)/cont.prop[,1]

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

      otuzero3.i <- (0.7754*otuzero.i)-4.2185
      otuzero3.i <- round(otuzero3.i)
      if(otuzero3.i <= 0){otuzero3.i <- 1}
      if(otuzero3.i > nrow(perc.diff.i)){otuzero3.i <- nrow(perc.diff.i)}
      if(otuzero3.i <= 0){otuzero3.i <- 1}
      otuzero3.i <- perc.diff.i[otuzero3.i,1]
      otuzero3.blank.i <- cont[cont$OTU == otuzero3.i,2]
      otuzero3.blank.ratios.i <- cont[,2]/otuzero3.blank.i
      otuzero3.correction.i <- cont[cont$OTU == otuzero3.i,3]
      otuzero3.correction.i <- otuzero3.blank.ratios.i*otuzero3.correction.i}else{otuzero3.correction.i <- cont.unlabeled[,2]}

      mean.i <- otuzero3.correction.i

  corrected <- cont[,3:ncol(cont)]-mean.i
  corrected[corrected < 0] <- 0
  corrected <- round(corrected)

  corrected <- cbind.data.frame(labels,cont[,2],corrected)

  colnames(corrected) <- colnames(cont)
  corrected <- rbind.data.frame(corrected,subset(contamination,blank==0))

  colnames(corrected) <- c("OTU",c(1:(ncol(corrected)-1)))
  
  corrected}