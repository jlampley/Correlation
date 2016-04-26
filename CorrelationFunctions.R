library(ape) #utility fns
library(geiger) #utilty fns
library(corHMM)
library(phylolm)
library(phytools)
library(rotl)

VisualizeData <- function(phy, data) {
  windows()
  plot (phy)
  print(str(phy))
  print(dim(data))
  windows()
  barplot(data)
  print(table(data))
  windows()
  hist(data[,1])
  windows()
  hist(data[,2])
  plot(phy) 
  tiplabels 
	#Important here is to LOOK at your data before running it. Any weird values? Does it all make sense? What about your tree? Polytomies?
}

CleanData <- function(phy, data) {
 
Nelumbotree <- treedata(phy, data)
return (Nelumbotree)
  
	#treedata() in Geiger is probably my favorite function in R.
}

RunContrasts <- function(phy, data, output.pdf="PIC.pdf") {
micpic <- pic(data[,1], phy, scaled=TRUE, var.contrasts=TRUE)
nicpic <- pic(data[,2], phy, scaled=TRUE, var.contrasts = TRUE)
  windows()
  plot(micpic[,1],sqrt(micpic[,2]))
  windows()
  plot(nicpic[,1],sqrt(nicpic[,2]))
  testmic <- cor.test(micpic[,1],sqrt(micpic[,2]))
  testnic <- cor.test(nicpic[,1],sqrt(nicpic[,2]))
  
  pic12 <- cbind(micpic[,1],nicpic[,1]) #positivizing above results
  pictable2 <- pic12
  for (i in seq(pic12[,1])){
    if (pic12[i,1]<0) {
      pictable2[i,] <- pic12[i,]*(-1)
    } 
  }  
  testposedpicmic <- cor.test(pictable2[,1],pictable2[,2])
    pdf (output.pdf)
    plot(pictable2[,1],pictable2[,2], xlab="PICtrait1", ylab="PICtrait2")
    abline(lm(pictable2[,1]~pictable2[,2]-1)) #-1 to force to the origin
    dev.off()
 # PSE? save the variances to standardize the contrasts
#test whether the stdzation was appropriate - consult the Garland paper
#no corr. between the standard dev. and the contrasts (take the sqrt of the variances
 #                                                     then positivize it)
#then do the corelation of the positivized things
#then plot(
 # then do Pagel.
#)
  	#Include here approaches to save plots, look at your data, regress through the origin, and return the results.
return(list(micpic=micpic, nicpic=nicpic, correlationtestm=testmic, 
            correlationtestn=testnic, testposedpicmic=testposedpicmic))
}

RunPagel94 <- function(phy, data) {
  Pagelmeth1 <- fitPagel(phy, data[,1], data[,2], method="ace")
  Pagelmeth2 <- fitPagel(phy, data[,1], data[,2], method="fitMK")
  Pagelmeth3 <- fitPagel(phy, data[,1], data[,2], method="fitDiscrete")
	#Calculate the rate estimates under each model, the likelihood of 
  #the data under each model, and the model averaged rates. 	Return the results
  return(list(Pagelmeth1=Pagelmeth1, Pagelmeth2=Pagelmeth2, Pagelmeth3=Pagelmeth3))
  }

Runphylolm <- function(phy, data) {
	data <- as.data.frame(data)
  linregrmodel <- phylolm(formula=data[,1]~data[,2], data=data, 
	 phy=phy, model="BM")
  return(linregrmodel)
}