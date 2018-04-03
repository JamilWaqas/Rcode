#Packages
library(MASS)
library(glmnet)
library(listdtr)
library(e1071)
library(svmpath)
################################################################################
#DATA processing
################################################################################
diabetes<-read.table("https://web.stanford.edu/~hastie/Papers/LARS/data64.txt",header = T)
diabetesY<-read.table("https://web.stanford.edu/~hastie/Papers/LARS/diabetes.data",header = T)
Y<-diabetesY[,11]
diabetes<-cbind(diabetes,Y)
smp_size <- floor(0.75*nrow(diabetes))
set.seed(123)
train_ind <- sample(seq_len(nrow(diabetes)), size = smp_size)
train <- diabetes[train_ind,]
test <- diabetes[-train_ind,]
mse <- function(error)
{
  (mean(error^2))
}

x.train<- as.matrix((train[,1:64]))
y.train<-as.matrix(train[,65])
x.test<-as.matrix(test[,1:64])
y.test<-as.matrix(test[,65])
x1<-x.train[,1:10]
x2<-x.test[,1:10]
x3<-cbind(x1[,2:4],x1[,7],x1[,9:10])
x4<-cbind(x2[,2:4],x2[,7],x2[,9:10])
##################################################################
# SLOG
################################################################
SLOG <- function(x,y,l,times,thresh,start=NULL){
  xtx<-crossprod(x)
  xty<-crossprod(x,y)
  p<-length(xty)
  n<-length(y)
  b.cur <- sign(xty)*l/p
  if(!is.null(start)) b.cur <- start
  b.old<-b.cur
  vin<-1:p
  temp<-rep(0,p)
  conv=FALSE
  k<-1
  while(conv==FALSE){
    b <- b.cur[vin]
    p <- length(b)
    B.inv <- diag(l/abs(b),nrow=p)
    b.new<-tcrossprod(chol2inv(chol(B.inv+xtx[vin,vin])),t(xty[vin]))
    temp[vin]<-as.vector(b.new)
    temp[abs(temp)<=thresh]<-0
    b.new<-temp
    vin<-which(b.new!=0)
    b.cur <- as.vector(b.new)
    conv<-(sqrt(sum((b.cur-b.old)^2))/sqrt(sum(b.old^2)))<times
    b.old<-b.cur
    k<-k+1
  }
  return(t(t(as.matrix(b.cur))))
}

#64 features
pr1<-SLOG(x.train,y.train,207,1e-10,1e-13)
re1<-x.test%*%pr1
mse(re1-y.test)#24728.15
mean(re1)#-3.515465
var(re1)#1349.959
#10 features
pr2<-SLOG(x1,y.train,609,1e-10,1e-16)
re2<-x2%*%pr2
mse(re2-y.test)#25745.58
mean(re2)# -1.603233
var(re2)#119.6921
#6 features
pr3<-SLOG(x3,y.train,71,1e-10,1e-16)
re3<-x4%*%pr3
mse(re3-y.test)#24626.66
mean(t(re3))#-4.176168
var(re3)#2383.072
#########################################################################################
#KSLOG
########################################################################################
RBFSLOG <- function(x,y,l,times,thresh,g,start=NULL){
  xtx<-radial.kernel(x,x,g)
  xty<-crossprod(x,y)
  n<-length(y)
  b.cur <- sign(y)*l/n
  if(!is.null(start)) b.cur <- start
  b.old<-b.cur
  vin<-1:n
  temp<-rep(0,n)
  conv=FALSE
  k<-1
  while(conv==FALSE){
    b <- b.cur[vin]
    t <- length(b)
    B.inv <- diag(l/abs(b),nrow=t)
    b.new<-ginv((B.inv+sqrt(B.inv)%*%xtx[vin,vin]%*%sqrt(B.inv)))
    b.new<-sqrt(B.inv)%*%b.new
    b.new<-b.new%*%sqrt(B.inv)
    b.new<-b.new%*%y[vin]
    temp[vin]<-as.vector(b.new)
    temp[abs(temp)<=thresh]<-0
    b.new<-temp
    vin<-which(b.new!=0)
    b.cur <- as.vector(b.new)
    conv<-(sqrt(sum((b.cur-b.old)^2))/sqrt(sum(b.old^2)))<times
    b.old<-b.cur
    k<-k+1
  }
  return(t(t(as.matrix(b.cur))))
}
#KSLOG 64 feature
pred1<-RBFSLOG(x.train,y.train,0.71,1e-10,1e-16,1.5)
k<-radial.kernel(x.train,x.test,1.5)
res1<-t(pred1)%*%k
mse(t(res1)-y.test)
mean(res1)
var(t(res1))

#KSLOG 10 features
pred2<-RBFSLOG(x1,y.train,0.71,1e-10,1e-16,21.7)
k1<-radial.kernel(x1,x2,21.7)
res2<-t(pred2)%*%k1
mse(res2-y.test)
plot(y.test)
points(t(res2),col="blue")
#MSE is 2716.402 RBF 10 features
mean(res2)
var(t(res2))

#KSLOG 6 features
#MSE is 2549.37 RBF 6 features
pred3<-RBFSLOG(x3,y.train,0.71,1e-10,1e-16,38.9)
k3<-radial.kernel(x3,x4,38.9)
res3<-t(pred3)%*%k3
mse(t(res3)-y.test)
mean(t(res3))
var(t(res3))
############################################################
#LASSO
###########################################################
#LASSO 64 features
cv.out = cv.glmnet(x.train, y.train, alpha = 1)
largelam = cv.out$lambda.1se
lasso.mod = glmnet(x.train, y.train, alpha = 1, lambda = largelam)
lasso.pred = predict(lasso.mod, s=largelam, newx = x.test)
cat(paste0("The test MSE for the lasso model using CV is: ",round(mean((y.test - lasso.pred)^2),2)))
#MSE is 3048.98

# LASSO 10 features
cv.out = cv.glmnet(x1, y.train, alpha = 1)
largelam = cv.out$lambda.1se
lasso.mod = glmnet(x1, y.train, alpha = 1, lambda = largelam)
lasso.pred = predict(lasso.mod, s=largelam, newx = x2)
cat(paste0("The test MSE for the lasso model using CV is: ",round(mean((y.test - lasso.pred)^2),2)))
#MSE is 3001.62

#LASSO with 6 features
cv.out = cv.glmnet(x3, y.train, alpha = 1)
largelam = cv.out$lambda.1se
lasso.mod = glmnet(x3, y.train, alpha = 1, lambda = largelam)
lasso.pred = predict(lasso.mod, s=largelam, newx = x4)
cat(paste0("The test MSE for the lasso model using CV is: ",round(mean((y.test - lasso.pred)^2),2)))
#MSE is 3032.44
mean(lasso.pred)
var(lasso.pred)
###################################################################
#SVM
###################################################################
#64 features
newdata<-data.frame(x.test)
X<-data.frame(cbind(y.train,x.train))
attach(X)
# traditional way: 64 features
tuneResult <- tune(svm, y.train~age+sex+bmi+map+tc+ldl+hdl+tch+ltg+glu+age.2+bmi.2+map.2+tc.2+ldl.2+hdl.2+tch.2+ltg.2+glu.2+
                     age.sex+ age.bmi+age.map+ age.tc +age.ldl+age.hdl+age.tch+age.ltg+age.glu+sex.bmi+sex.map+sex.tc+sex.ldl+
                     sex.hdl+sex.tch+sex.ltg+sex.glu+bmi.map+bmi.tc+bmi.ldl+bmi.hdl+bmi.tch+bmi.ltg+bmi.glu+
                     map.tc+map.ldl+map.hdl+map.tch+map.ltg+map.glu+tc.ldl+tc.hdl+tc.tch+tc.ltg+tc.glu+
                     ldl.hdl+ldl.tch+ldl.ltg+ldl.glu+hdl.tch+hdl.ltg+hdl.glu+tch.ltg+tch.glu+ltg.glu, data= X,kernel = "radial",
                   ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)),gamma=c(1.664592e-02,3.616560,5.459025,7.712658,9.264875,7.061480e-02,9.127456e-02,
                                                                                 7.450840e-02,3.589205,1.232052e-01,6.179574,7.573259e-03,3.009287e-01,20.8,17.5))

# test:
m<-tuneResult$best.model
pred<-predict (m, newdata)
svrMSE <- mse(t(t(pred))-as.matrix(y.test))
#MSE is 2950.32
mean(pred)
var(as.matrix(pred))

#10 features
tuneResult <- tune(svm, y.train~age+sex+bmi+map+tc+ldl+hdl+tch+ltg+glu, data= X,kernel = "radial",
                   ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)),gamma=c(1.664592e-02,3.616560e-04,5.459025e-03,7.712658e-06,9.264875e-01,7.061480e-02,9.127456e-02))
m<-tuneResult$best.model
pred<-predict (m, newdata)
svrMSE <- mse(t(t(pred))-as.matrix(y.test))
#MSE is 2950.32
mean(pred)
var(pred)

#6 features 
tuneResult <- tune(svm, y.train~sex+bmi+map+hdl+ltg+glu, data= X,kernel = "radial",
                   ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)),gamma=c(1.664592e-02,3.616560e-04,5.459025e-03,7.712658e-06,9.264875e-01,7.061480e-02,9.127456e-02))
m<-tuneResult$best.model
pred<-predict (m, newdata)
svrMSE <- mse(t(t(pred))-as.matrix(y.test))
#MSE is 2950.32
mean(pred)
var(pred)

###############################################################
#KRR
##############################################################
x<- x.train
y<-as.numeric(y.train)
#64 features
obj<-krr(x,y)
krrpred<-predict(obj,x.test)
error2<- krrpred - y.test
mse(error2)
mean(krrpred)
var(krrpred)
#10 features
#MSE is 2944.569 for 10 features
x<- x1
x.test<-x2
obj<-krr(x,y)
krrpred<-predict(obj,x.test)
error2<- krrpred - y.test
mean(krrpred)
var(krrpred)

#6 features
x<- x3
x.test<-x4
obj<-krr(x,y)
krrpred<-predict(obj,x.test)
error2<- krrpred - y.test
mean(krrpred)
var(krrpred)
########################################
#Following link has the code for HSIC-LASSO
#http://www.makotoyamada-ml.com/hsiclasso.html
#############################################



