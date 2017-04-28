#Data Load
####################################################################################################################################################
iris<-read.csv("iris.txt", header=F)
attach(iris)
parkinson<-read.csv("https://archive.ics.uci.edu/ml/machine-learning-databases/parkinsons/parkinsons.data")
attach(parkinson)
str(parkinson)
require(MASS)
#####################################################################################
#iris and parkinson, for one attribute
#####################################################################################
#for iris p=1
i1.iris<-as.data.frame(iris[,c(4,5)])
iris1.train<-i1.iris[1:70,]
i1test.1<-i1.iris[71:100,]
iris1.test<-i1test.1[,c(1)]
#####################################################################################
###############################################################
LDA1.iris<-function(train, test){
	
n <- dim(train)[1] 

p <- dim(train)[2]

i1labeli<-unique(i1.train[,2])
i1lab.set.i1<-which(i1.train[,2]==i1labeli[1])
i1lab.set.i2<-which(i1.train[,2]==i1labeli[2])
i1obs.set.i1<-i1.train[i1lab.set.i1,1]
i1obs.set.i2<-i1.train[i1lab.set.i2,1]	

labels<-i1labeli
labelclass.1<-i1lab.set.i1
labelclass.2<-i1lab.set.i2
observed.1<-i1obs.set.i1
observed.2<-i1obs.set.i2

p1.hat <- (1/n) * length(labelclass.1)
p2.hat <- (1/n) * length(labelclass.2)
m1.hat <- sum(observed.1)/length(labelclass.1)
m2.hat <- sum(observed.2)/length(labelclass.2)
sigma.hat <- ((t(observed.1) - m1.hat)%*%t(t(observed.1) - m1.hat) + (t(observed.2) - m2.hat)%*%t(t(observed.2) - m2.hat)) / (n-2)
Prediction <- rep(0,30)
delta <- as.matrix(test)%*%solve(sigma.hat)%*%(m2.hat - m1.hat) + as.numeric(log(p2.hat/p1.hat) - 0.5*(m1.hat + m2.hat)%*%solve(sigma.hat)%*%(m2.hat - m1.hat))

Prediction[which(delta > -1)] <- 1
return(list(prediction=Prediction))
}
#####################################################################################
i1lda<-LDA1.iris(train=iris1.train,test=iris1.test)

#####################################################################################
#for parkinson p=1
p1lab<-as.data.frame(parkinson[,c(18)])
p1lab1<-p1lab[1:120,]
p1lab2<-p1lab[121:195,]
p1.new<-as.data.frame(parkinson[,c(24)])
parkinson1.test<-p1.new[121:195,]
parkinson1.train<-cbind(p1.new[1:120,],p1lab1)
###############################################################
p1lda<-LDA1.parkinson<-function(train, test, no.test.obs){
	
n<-dim(train)[1] 
p<-dim(train)[2]

labelp1<-unique(p1.train[,2])
p1lab.set.p1<-which(p1.train[,2]==labelp1[1])
p1lab.set.p2<-which(p1.train[,2]==labelp1[2])
p1obs.set.p1<-p1.train[p1lab.set.p1,1]
p1obs.set.p2<-p1.train[p1lab.set.p2,1]

labels<-labelp1
labelclass.1<-p1lab.set.p1
labelclass.2<-p1lab.set.p2
observed.1<-p1obs.set.p1
observed.2<-p1obs.set.p2	

p1.hat<-(1/n)*length(labelclass.1)
p2.hat <-(1/n)*length(labelclass.2)

m1.hat<-sum(observed.1)/length(labelclass.1)
m2.hat<-sum(observed.2)/length(labelclass.2)

sigma.hat<-((t(observed.1) - m1.hat)%*%t(t(observed.1) - m1.hat) + (t(observed.2) - m2.hat)%*%t(t(observed.2) - m2.hat)) / (n-2)
Prediction <- rep(0,75)
delta <- as.matrix(test)%*%solve(sigma.hat)%*%(m2.hat - m1.hat) + as.numeric(log(p2.hat/p1.hat) - 0.5*(m1.hat + m2.hat)%*%solve(sigma.hat)%*%(m2.hat - m1.hat))

Prediction[which(delta > 0)] <- 1
return(list(prediction=Prediction))
}
####################################################################################################################################################
p1lda<-LDA1.parkinson(train=parkinson1.train,test=parkinson1.test)
########################################################################
#
##################################################################
#iris
####################################################################################################################################################
iris.train<-iris[1:70,]
iristest.1<-iris[71:100,]
iris.test<-iristest.1[,c(1,2,3,4)]
####################################################################################################################################################
LDA.iris<-function(train, test){
	
n <-dim(train)[1] 
p <-dim(train)[2]


ilabeli<-unique(train[,5])
ilab.set.i1<-which(train[,5]==i1labeli[1])
ilab.set.i2<-which(train[,5]==i1labeli[2])
iobs.set.i1<-train[i1lab.set.i1,1:4]
iobs.set.i2<-train[i1lab.set.i2,1:4]	

labels<-ilabeli
labelclass.1<-ilab.set.i1
labelclass.2<-ilab.set.i2
observed.1<-iobs.set.i1
observed.2<-iobs.set.i2

p1.hat<-(1/n)*length(labelclass.1)
p2.hat<-(1/n)*length(labelclass.2)

m1.hat<-colSums(observed.1)/length(labelclass.1)
m2.hat<-colSums(observed.2)/length(labelclass.2)

sigma.hat<-((t(observed.1) - m1.hat)%*%t(t(observed.1) - m1.hat) + (t(observed.2) - m2.hat)%*%t(t(observed.2) - m2.hat)) / (n-2)
Prediction <- rep(0,30)

delta <- as.matrix(test)%*%solve(sigma.hat)%*%(m2.hat - m1.hat) + as.numeric(log(p2.hat/p1.hat) - 0.5*(m1.hat + m2.hat)%*%solve(sigma.hat)%*%(m2.hat - m1.hat))

Prediction[which(delta > -1)]<-1 

return(list(prediction=Prediction))
}    
##################################################################
ilda<-LDA.iris(train=iris.train,test=iris.test)
##################################################################
#parkinson multi-atributes
lab<-as.data.frame(parkinson[,c(18)])
lab1<-lab[1:120,]
p.label<-lab[121:195,]
p.new<-as.data.frame(parkinson[,c(2:17,19:22)])
p.test<-p.new[121:195,]
p.train<-cbind(p.new[1:120,],lab1)
#############################################################
LDA.parkinson<-function(train, test){
	
n <-dim(train)[1] 
p <-dim(train)[2]


labelp<-unique(p.train[,21])
lab.set.p1<-which(p.train[,21]==labelp[1])
lab.set.p2<-which(p.train[,21]==labelp[2])
obs.set.p1<-p.train[lab.set.p1,1:20]
obs.set.p2<-p.train[lab.set.p2,1:20]

labels <-labelp
labelclass.1<-lab.set.p1
labelclass.2<-lab.set.p2
observed.1<-obs.set.p1
observed.2<-obs.set.p2


p1.hat<-(1/n)*length(labelclass.1)
p2.hat<-(1/n)*length(labelclass.2)

m1.hat<-colSums(observed.1)/length(labelclass.1)
m2.hat<-colSums(observed.2)/length(labelclass.2)
sigma.hat<-((t(observed.1) - m1.hat)%*%t(t(observed.1) - m1.hat) + (t(observed.2) - m2.hat)%*%t(t(observed.2) - m2.hat)) / (n-2)

Prediction<-rep(0,nrow(test))

delta <- as.matrix(test)%*%solve(sigma.hat)%*%(m2.hat - m1.hat) + as.numeric(log(p2.hat/p1.hat) - 0.5*(m1.hat + m2.hat)%*%solve(sigma.hat)%*%(m2.hat - m1.hat))

Prediction[which(delta > 0)]<-1 

return(list(prediction=Prediction))
}   


##################################################################
plda<-LDA.parkinson(train=p.train,test=p.test)
