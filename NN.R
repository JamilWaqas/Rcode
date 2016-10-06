##########################################################################
#reading iris data in R
##########################################################################
iris<-read.csv("iris.txt", header=F)
attach(iris)
str(iris)
summary(iris)
##########################################################################
#normalisation
##########################################################################
#normalising V1--V4 by creating a function "normalise"
normalise<-function(x) {
	
	return( (x-min(x) ) /(max(x)-min(x)))
}
##########################################################################
#A new data frame, which is normalised, excluding V5.
iris_new<-as.data.frame(lapply(iris[,c(1,2,3,4)],normalise))
##########################################################################
#Splitting data(train and test data sets)
##########################################################################
iris_train<-iris_new[1:70,]

iris_test<-iris_new[71:100,]

labledV5<-factor(iris$V5,levels=c(-1,1),labels=c("Setosa","Versicolor"))

iris_train_target<-(labledV5[1:70])

iris_test_target<-(labledV5[71:100])

###################################################################################################################################################
#K-NNcode
##########################################################################################################################################################
#Initial setup(by creating a function that calculates distance, by using a euclid matric)
###################################################################################################################################################
setup.knn<- function(x, iris_train, train_labels, k)
{

  matrix.x <- matrix(as.numeric(x), nrow=nrow(iris_train), ncol=length(x), byrow=T)

  matrix.x <- (abs(as.matrix(iris_train)-matrix.x))^2

  d <- sqrt((rowSums(matrix.x)))

  d <- data.frame(dist=d,label=train_labels)

  d <- (d[order(d$dist),])

d <- d[1:k,]

  setup <- names(sort(-table(d$label)))[1]
return(setup)

}
###################################################################################################################################################
#K-NN Algorithm applied
###################################################################################################################################################
knn <- function(iris_train, iris_test, train_labels, k, train_sample=NULL)

{
  

  train_sample <- nrow(iris_train)

  sub_sample <- sample(1:nrow(iris_train), train_sample, replace=F)

  iris_train <- iris_train[sub_sample,]

  train_labels <- train_labels[sub_sample]
  

  results <- apply(iris_test, 1, function(x) setup.knn(x, iris_train, train_labels, k))

  return(results)

}
##########################################################################
#Case where k=1,3(iris data)
###########################################################################k=1

#k=1
outputiris1<-knn(iris_train, iris_test, iris_train_target, k=1, train_sample=NULL)

#k=3
outputiris3<-knn(iris_train, iris_test, iris_train_target, k=3, train_sample=NULL)

##########################################################################
#Prediction(100% prediction for both)
##########################################################################
sum(outputiris1==iris_test_target)/length(iris_test_target)*100
sum(outputiris3==iris_test_target)/length(iris_test_target)*100
##########################################################################


