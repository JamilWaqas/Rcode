#############################################################################################################################################
                                                 #COMENTS
###################################################################################
#The Simple Approach: Estimates the model on a peice of data, and then computes the 1-step forecast on the remaining data.

#The code can further be shortened by making time series part as functions and using just one time-series and one AA code for both data, since the task was of analysis rather then programming,thus, I didn't devote time in doing that(copy paste is quicker :-)). 

###############################################################################################################################################
                                               #INSTRCTIONS
################################################################################################################################################

#time-series code some time to run, so thats why it's been commented out. The csv files were created and the data was saved and used for AA. 

#There are two ways to run the code:

#1) Piece by piece 

#2) Running the entire script

#If you chose option 1) then please download the packages before running the code, in option 2) packages should be installed by themself.

#Code is not extremly clean, there are bits which can be shortened, but at different points I've reffered the code to the parts of the project/dissertation. For example if at the left in block letters its says figure 1 then it mean that, the code will produce figure 1.

#################################################################################################################################################

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
#Average Behaviour Comparison FIGURE 13 CODE
###########################################################
dev.new()
#postscript("AAvsAvg.eps")
par(mfcol=c(1,2),oma = c(0, 0, 2, 0))
plot((AALMin - AvgLMin),main="Daily Minimum Temperature",ylim=c(-85000000,0),ylab="Loss",xlab="Days",type="l")
plot((AALMax - AvgLMax),ylim=c(-80000000,0),main="Daily Maximum Temperature",ylab="Loss",xlab="Days",type="l")
mtext("AA vs Average",outer = TRUE, cex = 1.5)
#dev.off()

######################################################################################################
#FIGURE 10 CODE
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
#FIGURE 11 CODE
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

#############################################
#FIGURE 7 CODE(The code is taken from an on-line resource and little changes are made)
############################################
#Model Selection plots

# the figures require ggplot2 library and
# all packages it depends on

# generate the x predictor
x <- runif(100,-2,2)
# generate the y response
y <- 2*x^3 + x^2 - 2*x +5 + rnorm(100)
xy <- data.frame(x=x, y=y)
# specify the maximum polynomial degree that will be explored
max.poly <- 7
 
# cretaing data.frame which will store model predictions
# that will be used for the smooth curves in Fig. 1
x.new <- seq(min(x), max(x), by=0.1)
degree <- rep(1:max.poly, each=length(x.new))
predicted <- numeric(length(x.new)*max.poly)
new.dat <- data.frame(x=rep(x.new, times=max.poly),
                      degree,
                      predicted)
 
# fitting lm() polynomials of increasing complexity
# (up to max.degree) and storing their predictions
# in the new.dat data.frame
for(i in 1:max.poly)
{
  sub.dat <- new.dat[new.dat$degree==i,]
  new.dat[new.dat$degree==i,3] <- predict(lm(y~poly(x, i)),
                                          newdata=data.frame(x=x.new))
}
 
# plotting the data and the fitted models

ovll<-ggplot()+geom_point(aes(x, y), xy, colour="darkgrey") + geom_line(aes(x, predicted, colour=as.character(degree)),new.dat)+scale_colour_discrete(name = "Degree")+labs(title="Fit")

# creating empty data.frame that will store
# AIC and BIC values of all of the models
AIC.BIC <- data.frame(criterion=c(rep("AIC",max.poly),
                                  rep("BIC",max.poly)),
                      value=numeric(max.poly*2),
                      degree=rep(1:max.poly, times=2))
 
# calculating AIC and BIC values of each model
for(i in 1:max.poly)
{
  AIC.BIC[i,2] <- AIC(lm(y~poly(x,i)))
  AIC.BIC[i+max.poly,2] <- BIC(lm(y~poly(x,i)))
}
 
# function that will perform the "leave one out"
# crossvalidation for a y~poly(x, degree) polynomial
crossvalidate <- function(x, y, degree)
{
  preds <- numeric(length(x))
    for(i in 1:length(x))
    {
        x.in <- x[-i]
        x.out <- x[i]
        y.in <- y[-i]
        y.out <- x[i]
        m <- lm(y.in ~ poly(x.in, degree=degree) )
        new <- data.frame(x.in = seq(-3, 3, by=0.1))
        preds[i]<- predict(m, newdata=data.frame(x.in=x.out))
    }
  # the squared error:
  return(sum((y-preds)^2))
}
 
# crossvalidating all of the polynomial models
# and storing their squared errors in
# the "a" object
a <- data.frame(cross=numeric(max.poly))
for(i in 1:max.poly)
{
  a[i,1] <- crossvalidate(x, y, degree=i)
}
 
# plotting AIC and BIC against model complexity
# (which is the polynomial degree)
AIC.plot <- qplot(degree, value, data=AIC.BIC,
                  geom="line", linetype=criterion) +
                  xlab("Polynomial degree") +
                  ylab("Criterion value") +
                  labs(title="Information theory & Bayes")+
                  geom_segment(aes(x=3, y=400,
                                   xend=3, yend=325),
                  arrow = arrow(length = unit(0.3, "cm"),
                 angle=20, type="closed")) +
                 theme(legend.position=c(0.8,0.5))


par(mfcol=c(1,2))


# plotting crossvalidated squared errors agains
# model complexity
cross.plot <- qplot(1:max.poly,cross, data=a, geom=c("line"))+
                    xlab("Polynomial degree") +
                    ylab("Squared error") +
                    geom_segment(aes(x=3, y=400,
                                     xend=3, yend=200),
                    arrow = arrow(length = unit(0.3, "cm"),
                    angle=20, type="closed")) +
                    labs(title="Crossvalidation")


dev.new()
plot_grid(ovll,AIC.plot,cross.plot, labels=c("","",""), ncol = 3, nrow = 1)

#######################################################################################################
#FIGURE 8 AND 9 CODE
######################################################################################################
gm<- seq(from=-1,to=1,by=0.1)
eeta<-0.5
gen1<-exp(-eeta*(1-gm)^2)
gen2<-exp(-eeta*(-1-gm)^2)
x<- c(4,0)
y<- c(0,4)
z<- c(1,2)
l<- seq(from=-0.7071068,to=0.7071068,by=0.1)

#######################################################################################################
dev.new()
plot((1-gm)^2,(-1-gm)^2,xlab=expression(paste("Loss when ", omega,"=1")),ylab=expression(paste("Loss when ", omega,"=-1")),type="l",xaxs="i", yaxs="i")
par(xaxs="i", yaxs="i")

dev.new()
plot(gen1,gen2,xlab="",ylab="",type="l",xaxs="i", yaxs="i")
par(xaxs="i", yaxs="i")

###############################################################
#FIGURE 5 AND 6 CODE
#############################################################
#Fourier Series And ARMA(2,2) plot, data and bit of code taken from stack-exchange then modified/tweaked according to my needs.

y <- c(11.622967, 12.006081, 11.760928, 12.246830, 12.052126, 12.346154, 12.039262, 12.362163, 12.009269, 11.260743, 10.950483, 10.522091,  9.346292,  7.014578,  6.981853,  7.197708,  7.035624,  6.785289, 7.134426,  8.338514,  8.723832, 10.276473, 10.602792, 11.031908, 11.364901, 11.687638, 11.947783, 12.228909, 11.918379, 12.343574, 12.046851, 12.316508, 12.147746, 12.136446, 11.744371,  8.317413, 8.790837, 10.139807,  7.019035,  7.541484,  7.199672,  9.090377,  7.532161,  8.156842,  9.329572, 9.991522, 10.036448, 10.797905)
t <- 18:65
ssp <- spectrum(y,plot=F)  
per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
reslm <- lm(y ~ sin(2*pi/per*t)+cos(2*pi/per*t))


rg <- diff(range(y))
ARMA<-auto.arima(y)
points<-as.numeric(fitted(ARMA))
dev.new()
plot(y~t,ylim=c(min(y)-0.1*rg,max(y)+0.1*rg),lty=c(8),main="Fourier series and ARMA(2,2)  approximation",xlab="x")
lines(fitted(reslm)~t,col=4,lty=2) 
reslm2 <- lm(y ~ sin(2*pi/per*t)+cos(2*pi/per*t)+sin(4*pi/per*t)+cos(4*pi/per*t))
lines(fitted(reslm2)~t,col=3)    
tt<-lm(y~points)
lines(fitted(tt)~t,col="red")
legend("bottomleft", legend = c("data points", "sin wave", "2 harmonic","ARMA(2,2)"), bty = "n",lwd = 2, cex = 1.2, col = c("black", "blue", "green","red"), lty = c(NA, 3, 1,1),pch=c(1,NA,NA,NA))
dev.new()
par(mfcol=c(2,1))
Acf(y,main="Acf of y")
Pacf(y,main="Pacf of y")
###################################################################################################
#FIGURE 3 CODE
###################################################################################################
#Auto-covariance theoretical vs sample
a<-0.6;N<-500;Nlags<-50
x<-rnorm(1)
X<-x
for(i in 1:N){
	x<-a*x+(rnorm(1)*sqrt(1-a^2))
X<-rbind(X,t(x))
}

x<-as.numeric(-51:50)

posautocov<-a^abs(0:50)
negautocov<-posautocov
autocov<-as.numeric(cbind(posautocov,negautocov))


sample_posautocov<-acf(X,Nlags,"covariance",plot=F)$acf
sample_negautocov<-sample_posautocov
sample_autocov<-as.numeric(cbind(sample_posautocov,sample_negautocov))

set1<-as.matrix(cbind(x,sample_autocov))
set2<-as.matrix(cbind(x,autocov))

subset1 <- subset(set1[,2], set1[,1] >= 0 & set1[,1] >= -50)
subset2 <- subset(set2[,2], set2[,1] >= 0 & set2[,1] >= -50)

dev.new()
plot(set1[52:102,1],subset1, ylim=c(-0.1,1),xlim=c(-50,50),type="l",xlab="lags",ylab="ACF",main="Sample vs theoretical ACF")
par(new=T)
plot(-set1[52:102,1],subset1,ylim=c(-0.1,1) ,xlim=c(-50,50),type="l",xlab="lags",ylab="ACF")
par(new=T)
plot(set2[52:102,1],subset2,ylim=c(-0.1,1), xlim=c(-50,50),type="l",xlab="lags",ylab="ACF",col="red")
par(new=T)
plot(-set2[52:102,1],subset2,ylim=c(-0.1,1) ,xlim=c(-50,50),type="l",xlab="lags",ylab="ACF",col="red")
legend("topright", legend = c("Sample ACF", "theoretical ACF"), bty = "n",lwd = 2, cex = 1.2, col = c("black","red"), lty = c(1,1))
####################################################################################################
#FIGURE 14 CODE
########################################################################################
dev.new()
par(mfrow=c(2,4))
Acf(max_temp_aus,main="Max temperature data")
Acf(min_temp_aus,main="Min temperature data")
Acf(diff(max_temp_aus),main="Differenced Max temperature data")
Acf(diff(min_temp_aus),main="Differenced Min temperature data")

Pacf(max_temp_aus,main="Max temperature data")
Pacf(min_temp_aus,main="Min temperature data")
Pacf(diff(max_temp_aus),main="Differenced Max temperature data")
Pacf(diff(min_temp_aus),main="Differenced Min temperature data")
###########################################################################
#EXAMPLE ON AA WORKINGS
###########################################################################
#AA Example Working
g0<-(1/3*exp(-2*(0.7^2)))+(1/3*exp(-2*(0.6^2)))+(1/3*exp(-2*(0.4^2)))
g1<-(1/3*exp(-2*(0.7-1)^2))+(1/3*exp(-2*(0.6-1)^2))+(1/3*exp(-2*(0.4-1)^2))
gam0<--0.5*log(g0)
gam1<--0.5*log(g1)
gam<-0.5-((gam1-gam0)/2)
w11<-1/3*exp(-2*(0.7-1)^2)
w12<-1/3*exp(-2*(0.6-1)^2)
w13<-1/3*exp(-2*(0.4-1)^2)

total<-w11+w12+w13

p21<-w11/total
p22<-w12/total
p23<-w13/total

###############################################################################
#FOLLOWING IS CODE FOR MATHEMATICA FIGURE 4 ( 14 DAYS TRIAL VERSION  USED)
###############################################################################
#AreaUnderTheCurve[f_, a_, b_, n_] := 
 #Module[{h = (b - a)/n, rects}, 
  #rects = Table[
   # Rectangle[{i, 0.}, {i + h, f[i + h/2]}], {i, a, b - h, h}];
 # Plot[f[x], {x, a, b}, 
  #  Epilog -> {EdgeForm[Black], FaceForm[None], rects}, 
   # ImageSize -> 500, Axes -> False, 
    #Frame -> {{True, False}, {False, True}}, 
    #FrameLabel -> {{"I"[\[Omega]], None}, {None, None}}, 
    #BaseStyle -> {FontSize -> 18}, FrameTicks -> None] Plot[
    #f[x], {x, a, b}, 
    #Epilog -> {EdgeForm[Black], FaceForm[None], rects, 
     # Arrow[{{a + 3 h, 2*f[a + 3 h]/5}, {a + 3.4 h, 2*f[a + 3 h]/5}}],
      # Arrow[{{a + 4 h, 2*f[a + 3 h]/5}, {a + 3.6 h, 
       #  2*f[a + 3 h]/5}}], 
      #Text[Style["\!\(\*FractionBox[\(2  \[Pi]\), \(T\)]\)", 
       # 14], {a + 3.5 h, 2 f[a + 3 h]/7}], 
      #Text[Style["0", 14], {0, -0.05}], 
     # Text[Style["\[Pi]", 14], {3, -0.05}], 
     # Text[Style["\!\(\*SubscriptBox[\(\[Omega]\), \(8\)]\)", 
      #  14], {a + 7.5 h, -0.05}], 
      #Text[Style["I[\!\(\*SubscriptBox[\(\[Omega]\), \(8\)]\)]", 
       # 14], {a + 7.5 h, f[a + 7.5 h] + 0.1}], Red, 
      #Line[{{a + 7.5 h, 0}, {a + 7.5 h, f[a + 7.5 h]}}]}, 
    #ImageSize -> 500, Axes -> False, 
    #Frame -> {{True, False}, {False, True}}, 
   # FrameLabel -> {{"I(\[Omega])", True}, {True, True}}, 
    #BaseStyle -> {FontSize -> 18}, FrameTicks -> None, 
    #PlotRangePadding -> {{Scaled[0.03], Scaled[0.03]}, {Scaled[0.07], 
     #  Scaled[0.03]}}]
  #]
#AreaUnderTheCurve[(2 #)/(1 + #^2) &, 0, Pi, 10]

#Code help by various personals on-line blogs/cites

###################################################################

#Remark:- In this developing phase, I was glad to understand mathematics and programming and tweak according to my needs. I look forward to understanding more complicated material in both deciplnes.
#####################################################################

BBMin<-AALossMinL
BBMax<-AALossMaxL

