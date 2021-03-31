rm(list=ls())
#library(nnls)
library(readr)
library(splines2)
library(splines)

library(foreach)
library(doParallel)

RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}
cleanData <- function(data , nbins){
  mydata <- data
  speedRange <- range(mydata$Speed)
  bins <- seq(from = speedRange[1], to = speedRange[2], by = ((speedRange[2] - speedRange[1])/(nbins- 1)) )
  binmeans <- seq(from = (bins[1]+bins[2])/2, to = (bins[nbins]+bins[nbins-1])/2, by = ((speedRange[2] - speedRange[1])/(nbins- 1)) ) 
  meanBin <- c()
  sdBin <- c()
  for ( i in 1:(nbins-1)){
    rangedata <- mydata[mydata$Speed <= bins[i+1],]
    rangedata <- rangedata[rangedata$Speed >= bins[i],]
    meanVal <- mean(rangedata$Power)
    meanBin <- c(meanBin , meanVal)
    sdVal <- sd(rangedata$Power)
    sdBin <- c(sdBin , sdVal)
  }
  se.bands=cbind(meanBin +3* sdBin ,meanBin -3* sdBin)
  
  removedData <- c()
  newdata <- c()
  for ( i in 1:(nbins-1)){
    
    rangedata <- mydata[mydata$Speed <= bins[i+1],]
    rangedata <- rangedata[rangedata$Speed >= bins[i],]
    n<- 1#dim(rangedata)[1]
    
    tempdata <- rangedata[(rangedata$Power - ( meanBin[i]- (3*sdBin[i]/(n^0.5)) ) ) >=0,]
    tempdata2 <- tempdata[(tempdata$Power - ( meanBin[i]+ (3*sdBin[i]/(n^0.5)) ) ) <=0,]
    
    newdata <- rbind(newdata , tempdata2)
    
    
    tempdata3 <- rangedata[(rangedata$Power - ( meanBin[i]- (3*sdBin[i]/(n^0.5)) ) ) < 0,]
    tempdata4 <- rangedata[(rangedata$Power - ( meanBin[i]+ (3*sdBin[i]/(n^0.5)) ) ) >0,]
    
    removedData <- rbind(removedData , tempdata4)
    removedData <- rbind(removedData , tempdata3)
    
  }
  cleanData <- newdata
}
# Function that returns Root Mean Squared Error
rmse <- function(error)
{
  #sqrt(mean(error^2))
  mean(abs(error))
}

# Function that returns Mean Absolute Error
mae <- function(error)
{
  mean(abs(error))
}

# Function that returns Mean Absolute Error
nmpae <- function(error , maxFitted)
{
  mean(abs(error)/maxFitted)*100
}

naToMean <- function (y){
  for(i in 1:length(y)){
    if(is.na(y[i])){
      if(i == 1){
        y[1] <- 0
      }else{
        y[i]<-y[i-1]
      }
    }
  }
  return(y)
}

match.numeric <- function(x, table) {
  are.equal <- function(x, y) isTRUE(all.equal(x, y))
  match.one <- function(x, table)
    match(TRUE, vapply(table, are.equal, logical(1L), x = x))
  vapply(x, match.one, integer(1L), table)
}
#######################
####Theoretical Data###
#######################

p<-c( 0, 0 ,  0, 0, 0, 100,  200,  400,  680,  950, 1200, 1400, 1535, 1610, 1650, 1650, 1650, 
      1650, 1650, 1650, 1650, 1650)
v<-c(0, 1, 2, 3, 3.5,  4.5,  5.0,  6.0,  7.0,  8.0 , 9.0, 10.0, 11.0, 12.0 ,13.0,
     14, 15, 16, 17, 18, 19, 20)

data.T.21<-data.frame(v, p)
v<-data.T.21$v
p<-data.T.21$p
uc<-3.5
ur<-13
us<-20
Speed.grid<-seq(0, 20, by=1)
fit=lm(p~poly(v,20),data=data.T.21)
theory_preds=predict(fit,newdata=list(v=Speed.grid),se=TRUE, rm.na=T)
########################
########################

#######################
####Config Variables###
#######################
showPlots = FALSE
tenMin = TRUE
outlierClean = FALSE
set.seed(127)
df = 7
cvNumbers<-5
tur = 33 #21, 49, 54 
coefWeigh =  1
coefWeighList = seq(0.0, 1, by=0.1)
bins_number = 10
######################
######################
target.curve = 2 #1 = Theoretical  # 2 = Farm AVG
if(target.curve==2){
  load( file = "C:/Statistic Mac/R code/Farm Average/FarmAverage.RData")
}
############################
####Load Data and Poolish###
############################
if(tenMin){
  mydata <- read_csv(paste("C:/Statistic Mac/R code/turbine data/Turbine/T",tur,".csv" ))
  mydata <- mydata[complete.cases(mydata), ]
  
}else{
  load(file = paste("/Users/mehrjoom/Desktop/R code/hybrid method/turbine data/hourly T",tur,".RData"))
  mydata <-mydata[complete.cases(mydata$Speed), ]
  mydata <-mydata[complete.cases(mydata$Power), ]
}
mydata <- mydata[mydata$Power >= 0, ]
mydata <- mydata[mydata$Power < 2000, ]
mydata <- mydata[mydata$Speed > 0, ]
mydata <- mydata[mydata$Speed <= 25, ]


mydata<-cleanData(mydata,10)
mydata$Speed = round(mydata$Speed*10)
mydata$Speed = mydata$Speed /10
if(outlierClean){
  y<-which(mydata$Speed > 12)
  u<-which(mydata$Power<1000)
  inds <-intersect(y,u)
  if(length(inds)>0){
    mydata<- mydata[-inds,]
  }
  y<-which(mydata$Speed < 5)
  u<-which(mydata$Power>1000)
  inds <-intersect(y,u)
  if(length(inds)>0){
    mydata<- mydata[-inds,]
  }
}
mydata <-mydata[complete.cases(mydata$Power), ]
mydata <-mydata[complete.cases(mydata$Speed), ]
newdata <- mydata
xlims =range(newdata$Speed)
maxPowerValue <- max(newdata$Power)

#x.grid=seq (from=xlims [1], to=xlims [2] , by=0.1)#length.out = ((xlims [2]-xlims [1])/0.1))
#xNumber = length(x.grid)
#new <- data.frame(x = x.grid)
#n <- length(x.grid)

x.grid=seq (from=0, to=20 , by=0.1)#length.out = ((xlims [2]-xlims [1])/0.1))
xNumber = length(x.grid)
new <- data.frame(x = x.grid)
n <- length(x.grid)


######################
######################

###########################
###Fit Natural Spline######
###########################

x<- newdata$Speed
fitNS <- lm(newdata$Power~ns(x ,df =df) ,data=newdata)
pred=predict (fitNS ,data.frame(x = new),se=T)
pred$fit[which(pred$fit<0)]=0
if(showPlots){
  plot(newdata$Speed, newdata$Power, pch=20 , col="grey")  
  lines(x.grid, pred$fit)
  if(target.curve == 2){
    lines(x.grid, avg.pred[0:length(x.grid)] , col = "green")
  }
}
##########################
##########################
##########################


##########################
###Calculate Error^2######
##########################
errors <- c()
varData <- c()
for( t in 1:dim(newdata)[1]){
  ind <- which(abs(x.grid - newdata[t,]$Speed)<0.0001)
  predictedVal <- pred$fit[ind]
  err <-  abs(newdata[t,]$Power - predictedVal)^2
  errors <- c(errors ,err)
  varData <- rbind(varData,c(newdata[t,]$Speed,err))
}
if(showPlots){
  plot(newdata$Speed, errors, pch=20 , col="grey")  
}
###Fit log(errors)~log(train$Speed)

##########################
###Calculate Weights######
##########################
if(length(which(errors==0))>0){
  error.speed = newdata$Speed[-which(errors==0)]
  errors = errors[-which(errors==0)]
}else{
  error.speed = newdata$Speed
}

fitVar <- lm(log(errors)~log(error.speed) ,data=as.data.frame(varData) )
beta <- fitVar$coefficients[2]
intercept <- fitVar$coefficients[1]
if(showPlots){
  plot(log(error.speed), log(errors), pch=20 , col="grey")  
  lines(log(x.grid), log(x.grid)*beta+intercept, col= "blue")
}
allW =  exp(1)^ (log(x.grid)*beta)

#Wnew= (1/allW) / sum(1/allW)#(max(allW))

if(showPlots){
  plot(x.grid, allW,xlab="Wind Speed", ylab="Weight" , type = "l")
}
###########################
###########################
###########################


########################################
####### Weight based on 1/SD ###########
########################################
bins = seq(from=xlims [1], to=xlims [2] , length.out = bins_number)#seq(from=xlims [1], to=xlims [2] , length.out = 10)
SDbins <- c()
x.grid.SD = rep(1,   length(x.grid))
i=1
for (bin in bins){
  datainBin <- newdata[which(bins[i]<=newdata$Speed & newdata$Speed<bins[i+1]),]
  sd1 <- sd(datainBin$Power)
  SDbins <- c(SDbins, sd1) 
  
  
  x.grid.SD[which(bins[i]<=x.grid & x.grid<bins[i+1])] = SDbins[i]
  i=i+1
  
}

wts = rep(1,   length(newdata$Speed))
bins = bins[complete.cases(SDbins)]
SDbins = SDbins[complete.cases(SDbins)]
target = SDbins
target[which(bins<10)] = min(target)
weight_data = c(bins)
weight_data =as.data.frame(t(rbind(weight_data,target)))
colnames(weight_data) <- c( "x","weight")
fitWeight <- lm(weight~ns(x ,df =(5)) ,data=weight_data)
weight_fitted=predict (fitWeight ,new,se=T)
for (sp in x.grid) {
  wts[which(abs(sp - newdata$Speed)<.0001)] = weight_fitted$fit[which(new$x==sp)]
}
if(showPlots){
  plot(x.grid , 1/x.grid.SD)
  #plot(x.grid , SDbins[which(bins)])
  plot(bins, target, col="black",  pch=20, 
       cex=0.5, xlab="Wind Speed (m/s)", ylab=" Weight",  
       xlim=c(0, 21),
       family="Times",
       cex.lab=1.25, 
       cex.axis=1)
  lines(x.grid , weight_fitted$fit)
  
  
}
weight_scheme2 = weight_fitted$fit + (0-min(weight_fitted$fit))
wiegh.x =  target + (0-min(weight_fitted$fit))
wiegh.x = wiegh.x / max(wiegh.x)
weight_scheme2 = weight_scheme2 / max(weight_scheme2)
if(showPlots){
  #plot(x.grid , SDbins[which(bins)])
  plot(bins, wiegh.x, col="black",  pch=20, 
       cex=0.5, xlab="Wind Speed (m/s)", ylab=" Weight",  
       xlim=c(0, 23),ylim = c(0,1),
       family="Times",
       cex.lab=1.25, 
       cex.axis=1)
  lines(x.grid, weight_scheme2)
}
########################
best_poly = c()
best_NSlist<- c()
best_NShyblist <- c()
best_SDNShyblist <- c()
best_lR_list = c()
best_lR_Hyb_list1 = c()
best_lR_Hyb_list2 = c()
df_list = seq(from = 4, to = 18 , by =2)
span_coef = 0.0125

ind_list = seq(from = 1 , to = nrow(newdata))
for (df  in df_list) {
  
  poly = c()
  NSlist<- c()
  NShyblist <- c()
  SDNShyblist <- c()
  LRlist = c()
  LRHyblist1 = c()
  LRHyblist2=c()
  targetList=c()
  for(coefWeigh in coefWeighList){
    HybridWeights = (coefWeigh*allW) / (max(allW))
    HybridWeights2 = coefWeigh*weight_scheme2 
    
    
    #######################
    ####setup   parallel###
    #######################
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    ########################
    ########################
    finalMatrix <- foreach(iterat=1:cvNumbers, .combine=cbind , .packages='splines') %dopar% {
     #for(iterat in 1:cvNumbers) {
      splits <-sample(nrow(newdata), replace=T)
      train <- newdata[unlist(splits),]#newdata[- unlist(splits[iterat]),]
      t=unlist(unique(splits))
      test <- newdata[ind_list[-t],]#newdata[unlist(splits[iterat]),]
      x<-train$Speed
      y<-train$Power
      data<-train
      
      ##### theori new version
      
      
      if(target.curve==1){
        theory_preds2=predict(fit,newdata=list(v=new$x),se=TRUE, rm.na=T)
        theory_preds2$fit[which(new$x > ur)] = maxPowerValue
        theory_preds2$fit[which(new$x < uc)] = 0
        target.fit = theory_preds2$fit
      }else if(target.curve == 2){
        target.fit = avg.pred
      }
      #######################
      #######################
      
      #######################
      ###Plynomial###########
      #######################
      fitpoly=lm(y~poly(x,df),data=data)
      predpoly=predict(fitpoly ,new,se=T)
      #######################
      ###Natural Spline######
      #######################
      fitNS <- lm(y~ns(x ,df =df) ,data=data)
      pred=predict (fitNS ,new,se=T)
      
      #######################
      ###Natural Spline######
      #######################
      fitLR <- loess(y~x ,data=data , span=(span_coef*df), degree = 0)
      predLR=predict (fitLR ,new,se=T)
      predLR$fit=naToMean(predLR$fit)
      ##############################
      #### Hybrid Natural Spline####
      ##############################
      LRhybFit = (HybridWeights*target.fit)+((1-HybridWeights)*predLR$fit)
      ##############################
      #### Hybrid Natural Spline####
      ##############################
      LRhybSDFit = (HybridWeights2*target.fit)+((1-HybridWeights2)*predLR$fit)
      
      ##############################
      #### Hybrid Natural Spline####
      ##############################
      NShybFit = (HybridWeights*target.fit)+((1-HybridWeights)*pred$fit)
      
      #################################################
      #### Hybrid Natural Spline Standard Deviation####
      #################################################
      NShybFitSD = (HybridWeights2*target.fit)+((1-HybridWeights2)*pred$fit)
      #plot(test$Speed, test$Power,col="grey",  pch=20)
      #lines(x.grid, pred$fit , col = "green")
      #lines(x.grid, NShybFit , col = "blue")
      #lines(x.grid, NShybFitSD , col = "red")
      #lines(x.grid, target.fit , col = "brown")
      #lines(x.grid, LRhybSDFit , col = "blue")
      #lines(x.grid, LRhybFit , col = "orange")
      #lines(x.grid, predLR$fit , col = "black")
      #######################
      #######################
      predpoly$fit[which(predpoly$fit<0)]=0
      pred$fit[which(pred$fit<0)]=0
      NShybFit[which(NShybFit<0)]=0
      NShybFitSD[which(NShybFitSD<0)]=0
      
      predLR$fit[which(predLR$fit<0)]=0
      LRhybSDFit[which(LRhybSDFit<0)]=0
      LRhybFit[which(LRhybFit<0)]=0
      
      LRhybFit[which(LRhybFit>maxPowerValue)]=maxPowerValue
      LRhybSDFit[which(LRhybSDFit>maxPowerValue)]=maxPowerValue
      predLR$fit[which(predLR$fit>maxPowerValue)]=maxPowerValue
      #######################
      #######################
      plerr = c()
      NSerr <- c()
      NSHyberr <- c()
      SDHybrerr = c()
      LRErr = c()
      LRHybErr1 = c()
      LRHybErr2 = c()
      targerErr=c()
      for( t in 1:dim(test)[1]){
        ind <- which(abs(x.grid - test[t,]$Speed)<.0001)
        
        # err <- (test[t,]$Power - target.fit[ind])
        # targerErr = c(targerErr,err)
        err <- (test[t,]$Power - predpoly$fit[ind])
        plerr <- c(plerr , err)


        err <- (test[t,]$Power - pred$fit[ind])
        NSerr <- c(NSerr , err)

        err <- (test[t,]$Power - NShybFit[ind])
        NSHyberr <- c(NSHyberr , err)


        err <- (test[t,]$Power - NShybFitSD[ind])
        SDHybrerr <- c(SDHybrerr , err)

        err <- (test[t,]$Power - predLR$fit[ind])
        LRErr <- c(LRErr , err)

        err <- (test[t,]$Power - LRhybFit[ind])
        LRHybErr1 <- c(LRHybErr1 , err)

        err <- (test[t,]$Power - LRhybSDFit[ind])
        LRHybErr2 <- c(LRHybErr2 , err)
      }
      #temp_var_value = (length(test$Power) * var(test$Power))
      # res = c( 1 - (  (sum( plerr ^ 2)) /  temp_var_value )  ) #c(rmse(plerr))
      # res = rbind(res ,1 - (  (sum( NSerr ^ 2)) /  temp_var_value )  ) #rmse(NSerr))
      # res = rbind(res , 1 - (  (sum( NSerr ^ 2)) /  temp_var_value )  ) #rmse(NSHyberr))
      # res = rbind(res , 1 - (  (sum( NSHyberr ^ 2)) /  temp_var_value )  ) # rmse(SDHybrerr))
      # res = rbind(res , 1 - (  (sum( LRErr ^ 2)) /  temp_var_value )  ) #rmse(LRErr))
      # res = rbind(res , 1 - (  (sum( LRHybErr1 ^ 2)) /  temp_var_value )  ) #rmse(LRHybErr1))
      # res = rbind(res , 1 - (  (sum( LRHybErr2 ^ 2)) /   temp_var_value ) ) #rmse(LRHybErr2))
      res = c(rmse(plerr))
      res = rbind(res ,rmse(NSerr))
      res = rbind(res , rmse(NSHyberr))
      res = rbind(res ,  rmse(SDHybrerr))
      res = rbind(res , rmse(LRErr))
      res = rbind(res , rmse(LRHybErr1))
      res = rbind(res , rmse(LRHybErr2))
      res
      #targetList = c(targetList , mae(targerErr))
    }
    
    stopCluster(cl)
    poly = c(poly , mean(finalMatrix[1,]))
    NSlist<- c(NSlist,mean(finalMatrix[2,]))
    NShyblist = c(NShyblist,mean(finalMatrix[3,]))
    SDNShyblist = c(SDNShyblist,mean(finalMatrix[4,]))
    
    LRlist<- c(LRlist,mean(finalMatrix[5,]))
    LRHyblist1 = c(LRHyblist1,mean(finalMatrix[6,]))
    LRHyblist2 = c(LRHyblist2,mean(finalMatrix[7,]))
  }
  
  best_poly = c(best_poly, min(poly))
  best_NSlist<- c(best_NSlist, min(NSlist))
  best_NShyblist <- c(best_NShyblist, min(NShyblist))
  best_SDNShyblist <- c(best_SDNShyblist, min(SDNShyblist))
  
  
  best_lR_list <- c(best_lR_list, min(LRlist))
  best_lR_Hyb_list1 <- c(best_lR_Hyb_list1, min(LRHyblist1))
  best_lR_Hyb_list2 <- c(best_lR_Hyb_list2, min(LRHyblist2))
}
best_poly 
best_NSlist
best_NShyblist 
best_SDNShyblist

# max(best_poly) 
# max(best_NSlist)
# max(best_NShyblist) 
# max(best_SDNShyblist)
# 
# max(best_lR_list) 
# max(best_lR_Hyb_list1) 
# max(best_lR_Hyb_list2) 

min(best_poly)
min(best_NSlist)
min(best_NShyblist)
min(best_SDNShyblist)

min(best_lR_list)
min(best_lR_Hyb_list1)
min(best_lR_Hyb_list2)
