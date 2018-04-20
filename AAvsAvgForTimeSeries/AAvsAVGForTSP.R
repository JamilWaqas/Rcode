
#Packages
install.packages("MASS")
install.packages("forecast")
install.packages("McSpatial")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("cowplot")
#Loading:
require(cowplot)
require(reshape2)
require(ggplot2)
require(MASS)
require(forecast)
require(McSpatial)

#Some pieces taken from stackoverflow
##########################################################################
#MAXIMUM DAILY TEMPERATURES DATA LOADING
##########################################################################

DailyTempAusMax<- read.csv("MaxTempAus.csv",T)
max_temp_aus<-ts(DailyTempAusMax,frequency=365,start=1980)

#########################################################################
#MAXIMUM DAILY TEMPERATURE PREDICTIONS FROM THE MODEL USED
#########################################################################

#Uncoment the following lines if you wish to run the code. The file PredictionsBigDataMax.csv reffers to these results.

#Cross-Validation for optimum value of k

#ymax<-window(max_temp_aus,frquency=365,start=1980,end=c(1981,0))
#xmax<-1:length(ymax)
#qmax<-floor(365/2)-1
#fit_gmax<-McSpatial::fourier(ymax~xmax,minq=1,maxq=qmax-1,crit="gcv")
#fit_gmax$q

#Simple approach:
#p<-0:2;d<-0:1;q<-0:2
#combmax<-as.matrix(expand.grid(p,d,q))
#fcmax<-matrix(0,nrow=3285,ncol=nrow(comb))
#for (k in 1:(nrow(comb))){
 #    p<- comb[k,1];d<- comb[k,2];q<- comb[k,3]
 #    trainmax<-window(max_temp_aus,end=1980.999)
  #   fitmax<-Arima(trainmax, order=c(p,d,q),method="ML",xreg=forecast::fourier(trainmax,95))
  #   refitmax <-Arima(max_temp_aus, model=fitmax,xreg=forecast::fourier(max_temp_aus,95))
  #   fcmax[,k] <-window(fitted(refitmax),start=1981)
  #  }      
#write.csv(fcmax,"PredictionsBigDataMax.csv",row.names=F)

#######################################################################################
#IMPLEMENTATION OF AA ON DAILY MAXIMUM TEMPERATURE DATA
#########################################################################################
#Algorithm preperation:
pre_max<- t(read.csv("PredictionsBigDataMax.csv",T))
pre_aus_max<- DailyTempAusMax[366:3650,1]
expertsPredictionsMax<- pre_max
outcomesMax<- t(pre_aus_max)
N<-nrow(expertsPredictionsMax)
col<-ncol(expertsPredictionsMax)
row<- nrow(expertsPredictionsMax)
AMax<-min(outcomesMax)
BMax<-max(outcomesMax)
etaMax<- 2 / ((BMax-AMax)^2)

#Substitition Function:
substitutionFunction<-function(p,ep){
	gAMax<- -(1/etaMax) * log(p %*% exp(-(etaMax) * t(ep - AMax)^2))
	gBMax<- -(1/etaMax) * log(p %*% exp(-(etaMax) * t(ep - BMax)^2))
 gammaMax<- (0.5*(BMax + AMax)) - ((gBMax - gAMax)/(2 * (BMax - AMax)))
	return(gammaMax)
}

#Aggregation Algorithm:
AAgpredictions<-function(expertsPredictions,outcomes){
	weights<-matrix(1,N,1)
	AApredictions<-matrix(1,col,1)
	for(t in 1:col){
	normalisedWeights<-weights/sum(weights) 
	    AApredictions[t]<-substitutionFunction(t(normalisedWeights),t(expertsPredictions[,t]))
	    weights1<-(normalisedWeights) * as.vector(exp(-etaMax * (expertsPredictions[,t] - outcomes[,t])^2))	
		weights<-weights1/sum(weights)
				}
		return(AApredictions)
}
predMax<-AAgpredictions(expertsPredictionsMax,outcomesMax)


outMax<- t(outcomesMax)
expertMax<-t(expertsPredictionsMax)
ExpertLossMax<-matrix(0,nrow=3285,ncol=18)
for(i in 1:18){
ExpertLossMax[,i]<-cumsum((expertMax[,i]-outMax)^2)
}
AALossMaxL<-cumsum((predMax - outMax)^2)
AvgLossMax<-cumsum(((as.matrix(rowSums(expertMax)/18)) - outMax)^2)
AALMax<-as.matrix(AALossMax)
AvgLMax<-as.matrix(AvgLossMax)

TotalAALossMax<-AALMax[3285,]
TotalAALossMax
TotalAvgLossMax<-AvgLMax[3285,]
TotalAvgLossMax
AggregationMinusExpertLossMax<-AALMax[3285,]-ExpertLossMax[3285,]
AggregationMinusExpertLossMax



#########################################################################
#MINIMUM DAILY TEMPERATURE PREDICTIONS FROM THE MODEL USED
#########################################################################


DailyTempAusMin<- read.csv("TempAus.txt",T)
conv<-as.numeric(as.character(unlist(DailyTempAusMin[[1]])))

#Data:
min_temp_aus<-ts(conv,frequency=365,start=1980)

#Uncoment the following lines if you wish to run the code. The file PredictionsBigDataMin.csv reffers to these results.


#Cross-Validation for optimum value of k
#y<-window(min_temp_aus,frquency=365,start=1980,end=c(1981,0))
#x<-1:length(y)
#qmax<-floor(365/2)-1
#fit_g<-McSpatial::fourier(y~x,minq=1,maxq=qmax-1,crit="gcv")
#fit_g$q

#Simple approach:

#p<-0:2;d<-0:1;q<-0:2
#comb<-as.matrix(expand.grid(p,d,q))
#fc<-matrix(0,nrow=3285,ncol=nrow(comb))
#for (k in 1:(nrow(comb))){
 #   p<- comb[k,1];d<- comb[k,2];q<- comb[k,3]
 #   train<-window(min_temp_aus,end=1980.999)
  #   fit<-Arima(train, order=c(p,d,q),method="ML",xreg=forecast::fourier(train,57))
  #   refit <-Arima(min_temp_aus, model=fit,xreg=forecast::fourier(min_temp_aus,57))
  #   fc[,k] <-window(fitted(refit),start=1981)
  #  }      
#write.csv(fc,"PredictionsBigDataMin.csv",row.names=F)


#######################################################################################
#IMPLEMENTATION OF AA ON DAILY Mminimum TEMPERATURE DATA
#########################################################################################

#Algorithm preperation:
pre_min<- t(read.csv("PredictionsBigDataMin.csv",T))
pre_aus_min<- min_temp_aus[366:3650]#read.csv("outcomes.csv",T)
expertsPredictionsMin<- pre_min
outcomesMin<- t(pre_aus_min)
N<-nrow(expertsPredictionsMin)
col<-ncol(expertsPredictionsMin)
row<- nrow(expertsPredictionsMin)
AMin<-min(outcomesMin)
BMin<-max(outcomesMin)
etaMin<- 2 / ((BMin-AMin)^2)

#Substitition Function:
substitutionFunction<-function(p,ep){
	gAMin<- -(1/etaMin) * log(p %*% exp(-(etaMin) * t(ep - AMin)^2))
	gBMin<- -(1/etaMin) * log(p %*% exp(-(etaMin) * t(ep - BMin)^2))
 gammaMin<- (0.5*(BMin + AMin)) - ((gBMin - gAMin)/(2 * (BMin - AMin)))
	return(gammaMin)
}

#Aggregation Algorithm:
AAgpredictions<-function(expertsPredictions,outcomes){
	weights<-matrix(1,N,1)
	AApredictions<-matrix(1,col,1)
	for(t in 1:col){
	normalisedWeights<-weights/sum(weights) 
	    AApredictions[t]<-substitutionFunction(t(normalisedWeights),t(expertsPredictions[,t]))
	    weights1<-(normalisedWeights) * as.vector(exp(-etaMin * (expertsPredictions[,t] - outcomes[,t])^2))	
		weights<-weights1/sum(weights)
				}
		return(AApredictions)
}
predMin<-AAgpredictions(expertsPredictionsMin,outcomesMin)

#Plots:
outMin<- t(outcomesMin)
expertMin<-t(expertsPredictionsMin)
ExpertLossMin<-matrix(0,nrow=3285,ncol=18)
for(i in 1:18){
ExpertLossMin[,i]<-cumsum((expertMin[,i]-outMin)^2)
}
AALossMinL<-cumsum((predMin - outMin)^2)
AvgLossMin<-cumsum(((as.matrix(rowSums(expertMin)/18)) - outMin)^2)
AALMin<-as.matrix(AALossMin)
AvgLMin<-as.matrix(AvgLossMin)


#write.csv(exp2,"regret2.csv",row.names=F)

TotalAALossMin<-AALMin[3285,]
TotalAALossMin
TotalAvgLossMin<-AvgLMin[3285,]
TotalAvgLossMin
AggregationMinusExpertLossMin<-AALMin[3285,]-ExpertLossMin[3285,]
AggregationMinusExpertLossMin

###############################################################################################
#IMPLEMENTATION WITH A DIFFERENT APPROACH(stack-overflow help for the 2 implementation also)
###############################################################################################



#set.seed(1234)
#y <- ts(sort(rnorm(30)), start = 1978, frequency = 1) # annual data
#fcasts <- numeric(10)
#train <- window(y, end = 1997) 
#fit<-arima(train)
#models<-list()
#for (i in 1:10) { # start rolling forecast
  # start from 1997, every time one more year included
 # win.y <- window(y, end = 1997 + i) 
  #refit <- Arima(win.y,model=fit)
 # fcasts[i] <- forecast(fit, h = 1)$mean
#}
#train <- window(y,end=1997)
#fit <- arima(train)
#refit <- Arima(y, model=fit)
#fc <- window(fitted(refit), start=1998)

#fc-fcasts# Same till the 16th decimal place!
##########################################
#FIGURE 12 CODE(stack-overflow help for this plot)
##########################################
dev.new()
#postscript("AAvsOut.eps")
par(mfcol=c(1,2),oma = c(0, 0, 2, 0))
plot(outMin,main="Daily Minimum Temperature",ylab="Temperature",xlab="Days")
points(predMin,col="red")
plot(outMax,main="Daily Maximum Temperature",ylab="Temperature",xlab="Days")
points(predMax,col="red")
mtext("AA vs Outcomes",outer = TRUE, cex = 1.5)
#dev.off()
###########################################################
#Average Behaviour Comparison 
###########################################################
dev.new()
#postscript("AAvsAvg.eps")
par(mfcol=c(1,2),oma = c(0, 0, 2, 0))
plot((AALMin - AvgLMin),main="Daily Minimum Temperature",ylim=c(-85000000,0),ylab="Loss",xlab="Days",type="l")
plot((AALMax - AvgLMax),ylim=c(-80000000,0),main="Daily Maximum Temperature",ylab="Loss",xlab="Days",type="l")
mtext("AA vs Average",outer = TRUE, cex = 1.5)
#dev.off()

######################################################################################################

######################################################################################################
#ggplot2 Plotting here(Was a bit different then using ususal plot or matplot in R, but found that it is so cool)
#############################################################################################################################
exp1<-ExpertLossMin-rep(AALossMin,18)
exp2<-ExpertLossMax-rep(AALossMax,18)

av1<-ExpertLossMin-rep(AvgLMin,18)
av2<-ExpertLossMin-rep(AvgLMax,18)
#install.packages("cowplot")
library(cowplot)

Days<-1:3285


#Min Temperature

df<-data.frame(cbind(Days,exp1[,1],exp1[,2],exp1[,3],exp1[,4],exp1[,5],exp1[,6],exp1[,7],exp1[,8],exp1[,9],exp1[,10],exp1[,11],exp1[,12],exp1[,13],exp1[,13],exp1[,15],exp1[,16],exp1[,17],exp1[,18]))
colnames(df)<-c("Days","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15","E16","E17","E18")


df.long<-melt(df,id.vars="Days")
colnames(df.long)<-c("Days","Expert","prediction")


sp<-ggplot(df.long,aes(Days,prediction,color=Expert))+geom_line()+labs(x="Days",y="Loss")+theme(legend.position='none') 

bp<-ggplot(df.long,aes(Days,prediction,color=Expert))+geom_line()+coord_cartesian(ylim=c(-2000, 4000))+labs(x="Days",y="Loss")+theme(legend.position='bottom')
fp<-bp + geom_hline(yintercept=-1061.359, linetype="dashed", color = "red")


#Max Temperature

df1<-data.frame(cbind(Days,exp2[,1],exp2[,2],exp2[,3],exp2[,4],exp2[,5],exp2[,6],exp2[,7],exp2[,8],exp2[,9],exp2[,10],exp2[,11],exp2[,12],exp2[,13],exp2[,13],exp2[,15],exp2[,16],exp2[,17],exp2[,18]))
colnames(df1)<-c("Days","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15","E16","E17","E18")
df1.long<-melt(df1,id.vars="Days")
colnames(df1.long)<-c("Days","Expert","prediction")

sp1<-ggplot(df1.long,aes(Days,prediction,color=Expert))+geom_line()+labs(x="Days",y="Loss")+theme(legend.position='none')
bp1<-ggplot(df1.long,aes(Days,prediction,color=Expert))+geom_line()+coord_cartesian(ylim=c(-2000, 4000))+labs(x="Days",y="Loss")+theme(legend.position='bottom')
fp1<-bp1 + geom_hline(yintercept=-1904.307, linetype="dashed", color = "red")
dev.new()
#setEPS()
#postscript("LossPlots.eps")
plot_grid(fp, fp1,sp, sp1, labels=c("A","B","C","D"), ncol = 2, nrow = 2)
#dev.off()
##################################################################

###################################################################
#Average 
df11<-data.frame(cbind(Days,av1[,1],av1[,2],av1[,3],av1[,4],av1[,5],av1[,6],av1[,7],av1[,8],av1[,9],av1[,10],av1[,11],av1[,12],av1[,13],av1[,13],av1[,15],av1[,16],av1[,17],av1[,18]))
colnames(df11)<-c("Days","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15","E16","E17","E18")

df11.long<-melt(df11,id.vars="Days")
colnames(df11.long)<-c("Days","Expert","prediction")


sp11<-ggplot(df11.long,aes(Days,prediction,color=Expert))+geom_line()+labs(x="Days",y="Loss")+theme(legend.position='none') 

bp11<-ggplot(df11.long,aes(Days,prediction,color=Expert))+geom_line()+coord_cartesian(ylim=c(-1000, -150000))+labs(x="Days",y="Loss")+theme(legend.position='bottom')
fp11<-bp11 + geom_hline(yintercept=-1061.359, linetype="dashed", color = "red")
###########################
#Max Temperature
##########################
df22<-data.frame(cbind(Days,av2[,1],av2[,2],av2[,3],av2[,4],av2[,5],av2[,6],av2[,7],av2[,8],av2[,9],av2[,10],av2[,11],av2[,12],av2[,13],av2[,13],av2[,15],av2[,16],av2[,17],av2[,18]))
colnames(df22)<-c("Days","E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15","E16","E17","E18")
df22.long<-melt(df22,id.vars="Days")
colnames(df22.long)<-c("Days","Expert","prediction")

sp22<-ggplot(df22.long,aes(Days,prediction,color=Expert))+geom_line()+labs(x="Days",y="Loss")+theme(legend.position='none')
bp22<-ggplot(df22.long,aes(Days,prediction,color=Expert))+geom_line()+coord_cartesian(ylim=c(-1000, -150000))+labs(x="Days",y="Loss")+theme(legend.position='bottom')
fp22<-bp22 + geom_hline(yintercept=-1904.307, linetype="dashed", color = "red")

#setEPS()
#postscript("AvgPlots.eps")
dev.new()
plot_grid(fp11, fp22,sp11, sp22, labels=c("E","F","G","H"), ncol = 2, nrow = 2)
#dev.off()
