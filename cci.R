#x<-c(17,20,10,17,12,15,19,22,17,19,14,22,18,17,13,12,18,15,17)
#Independence condition
set.seed(1)
cci<- function(x){
  ss<- sum(x)
  n<- length(x)+1
  mi<- min(x)
  ma<- max(x)
  if(ss > 0 | (ss == 0 && -mi <= -ma)){
    lower<- mi
    upper<- - (n*mi-(2*ss)) / (n-2)
  }else if(ss < 0 | (ss == 0 && -ma <= -mi) ){
    lower<- - (n*ma-(2*ss)) / (n-2)
    upper<- ma
  }
  return(list(lower_cci=lower,upper_cci=upper))
}
#violating independence/normality assumption
au<-numeric()
al<-numeric()
bu<-numeric()
bl<-numeric()
c<-numeric()
test<- numeric()
for(i in 1:1500){
  x<-arima.sim(list(order=c(1,0,0), ar=.5), n=25)
  x<-sample(x,20,F)
  test[i]<-shapiro.test(x)$p.value
  au[i]<-cci(x[1:19])$upper
  al[i]<-cci(x[1:19])$lower
  bu[i]<-t.test(x)$conf.int[2]
  bl[i]<-t.test(x)$conf.int[1]
  c[i]<-x[20]
}

CCI<-cbind(as.matrix(au),as.matrix(al),as.matrix(c),as.matrix(test))
tab_cci<-ifelse(CCI[,3]>CCI[,1]|CCI[,3]<CCI[,2],1,0)
table(tab_cci)

CI<-cbind(as.matrix(bu),as.matrix(bl),as.matrix(c),as.matrix(test))
tab_ci<-ifelse(CI[,3]>CI[,1]|CI[,3]<CI[,2],1,0)
table(tab_ci)

normality<- ifelse(CI[,4]<0.05,1,0)
table(normality)


au<-numeric()
al<-numeric()
bu<-numeric()
bl<-numeric()
c<-numeric()
test<- numeric()
for(i in 1:1500){
  x<-runif(25,0,1)
  x<-sample(x,20,F)
  test[i]<-shapiro.test(x)$p.value
  au[i]<-cci(x[1:19])$upper
  al[i]<-cci(x[1:19])$lower
  bu[i]<-t.test(x)$conf.int[2]
  bl[i]<-t.test(x)$conf.int[1]
  c[i]<-x[20]
}
CCI<-cbind(as.matrix(au),as.matrix(al),as.matrix(c),as.matrix(test))
tab_cci<-ifelse(CCI[,3]>CCI[,1]|CCI[,3]<CCI[,2],1,0)
table(tab_cci)

CI<-cbind(as.matrix(bu),as.matrix(bl),as.matrix(c),as.matrix(test))
tab_ci<-ifelse(CI[,3]>CI[,1]|CI[,3]<CI[,2],1,0)
table(tab_ci)

normality<- ifelse(CI[,4]<0.05,1,0)
table(normality)

au<-numeric()
al<-numeric()
bu<-numeric()
bl<-numeric()
c<-numeric()
test<- numeric()
for(i in 1:1500){
  x<-rnorm(25,0,1)
  x<-sample(x,20,F)
  test[i]<-shapiro.test(x)$p.value
  au[i]<-cci(x[1:19])$upper
  al[i]<-cci(x[1:19])$lower
  bu[i]<-t.test(x)$conf.int[2]
  bl[i]<-t.test(x)$conf.int[1]
  c[i]<-x[20]
}
CCI<-cbind(as.matrix(au),as.matrix(al),as.matrix(c),as.matrix(test))
tab_cci<-ifelse(CCI[,3]>CCI[,1]|CCI[,3]<CCI[,2],1,0)
table(tab_cci)

CI<-cbind(as.matrix(bu),as.matrix(bl),as.matrix(c),as.matrix(test))
tab_ci<-ifelse(CI[,3]>CI[,1]|CI[,3]<CI[,2],1,0)
table(tab_ci)

normality<- ifelse(CI[,4]<0.05,1,0)
table(normality)

