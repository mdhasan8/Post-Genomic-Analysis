setwd("C:/Users/Easin/OneDrive/Documents")
hcc<-read.csv("classify_HCC.csv")
data<-read.csv("gene_data.csv")
y<-factor(hcc$phenotype.tissue)
y<-as.numeric(y)
y<-ifelse(y==1,0,1)
y<-y[-c(435:445)]
v<-(data[1,])
dat<-data[-c(1,436:446),]
datf<-cbind(dat,y)
library(dplyr)
datfm <- sample_n(datf, 500)
y.m <- datfm $ y
datc <- sample(datfm,30)
colnames(datc)
datk <- cbind(datc, y.m)
colnames(datk)
dim(datk)

#partitioning the data#
n1<- NROW(datk)
set.seed(24)
id.test <- sample(1:2, size=n1, replace=TRUE,prob=c(2/3,1/3)) 
dat.training <- datk[id.test==1, ] 
dat.test<-datk[id.test==2,] 
dim(dat.training );dim(dat.test)
train<-dat.training[,-c(31)]
y.train<-dat.training[,c(31)]
test<-dat.test[,-c(31)]
y.test<-dat.test[,c(31)]

colnames(dat.training)
dat.training $ y.m
class(dat.training)
class(y.m)

miss.info<- function(dat, filename=NULL){ 
  vnames <- colnames(dat); vnames
  n <- nrow(dat) 
  out <- NULL 
  for (j in 1: ncol(dat)){ 
    vname <- colnames(dat)[j] 
    x <- as.vector(dat[,j])
    n1 <- sum(is.na(x), na.rm=T) 
    n2 <- sum(x=="NA", na.rm=T) 
    n3 <- sum(x=="", na.rm=T) 
    nmiss <- n1 + n2 + n3 
    ncomplete <- n-nmiss
    out <- rbind(out, c(col.number=j, vname=vname, mode=mode(x), n.levels=length(unique(x)), ncomplete=ncomplete, miss.perc=nmiss/n)) }
  out <- as.data.frame(out) 
  row.names(out) <- NULL 
  return(out)
} 

miss.info(dat.training)


dat.training[] <- lapply(dat.training, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(dat.training, class)

test[] <- lapply(test, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})
sapply(test, class)

library("e1071")
library("svmpath")
library("kernlab")
library(verification)
msvm<- ksvm(y.m ~ ., data = dat.training, 
            type = "C-bsvc", kernel = "rbfdot",
            kpar = list(sigma = 0.1), C = 10, 
            scaled=TRUE, prob.model = TRUE)

s.pred<-predict(msvm, test, type = "probabilities")
s.pred <- s.pred[, 1]
pred.svm<-ifelse(s.pred<=0.5,1,0)



#accuracy of SVM method 
accuracy.rf<- function(truth, predicted)
  if(length(truth) > 0)
    sum(truth==predicted)/length(truth) else    
      return(0)
accuracy.rf(y.test,hat.rf)

#sensitivity of SVM method
sensitivity.rf <- function(truth, predicted)
  # 1 means positive (present)
  if(sum(truth==1) > 0)
    sum(predicted[truth==1]==1)/sum(truth==1)   else
      return(0)
sensitivity.rf(y.test,hat.rf)

#specificity of SVM method
specificity.rf<- function(truth, predicted)
  if(sum(truth==0) > 0)
    sum(predicted[truth==0]==0)/sum(truth==0)   else
      return(0)
specificity.rf(y.test,hat.rf)

#area under roc curve
y.m.test <- y.test 

f.AUC.s<- roc.area(obs=y.test , pred=pred.svm)$A
f.L.s<- verify(obs=y.m.test, pred=pred.svm)
roc.plot(f.L.s, plot.thres = NULL, col="red")
text(x=0.7, y=0.2, paste("ROC =", round(f.AUC.s, digits=4), 
                         sep=" "), col="blue", cex=1.2)
f.AUC.s


# VARIABLE IMPORTANCE RANKING
round(importance(fit.rf), 2)
varImpPlot(fit.rf, main="Variable Importance Ranking")
# error rate based on training set
fit.rf$confusion
