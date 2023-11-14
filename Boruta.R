#Read the data
library(lessR)
library(tm)
library(topicmodels)
library(combinat)
setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/Thesis_NVP/Chapter1/Paper_Github/")
load(file = "Data.rda")
head(Data)

#class labels associated with the data
load(file = "actual.rda")
summary(actual)
BarChart(actual,ylab =" Count-Health Status",xlab="Health Status")

########################## Tau-Path Method ########################

tot_sub<-rowSums(Data)
#Calculating Abundance
Abnd_Data<-Data
for(i in 1:nrow(Data))
{
  Abnd_Data[i,]<-(Data[i,]/tot_sub[i])*100
}

Abnd_Measure<-colSums(Abnd_Data)
th<-1

#Calculating Prevelance
Prev_Data<-ifelse(Abnd_Data<=1,0,1)
Prev_Measure<-colSums(Prev_Data)


#Selecting top-K

cmat <- function(x,y){
  if(!is.numeric(x)) stop("x should be in numeric format \n")
  if(!is.numeric(y)) stop("y should be in numeric fomrat \n")
  
  n <- length(x)
  cmatrix <- matrix(0,nrow=n,ncol=n)
  
  for ( i in 1:n){
    for ( j in 1:n) {
      if((x[i]-x[j])*(y[i]-y[j]) > 0) cmatrix[i,j] <- 1
      else if((x[i]-x[j])*(y[i]-y[j]) < 0) cmatrix[i,j] <- -1
    }
  }
  
  return(cmatrix)
}

## A function to calculate tau-path ##
TauPath <- function(mat){
  ###input: C Matrix
  
  path <- NULL
  t <- numeric(0)
  for( i in 1:dim(mat)[1]-1){
    offave<- sum(mat[1:(1+i),1:(1+i)])/((1+i)^2-(1+i))
    t <- cbind(t,offave)
  }
  colnames(t) <- NULL
  t <- t[-1]
  return(t)
}

## A function to re-order our observations ##
ReOrder <- function(x,y){
  ##position of paris
  poi <- c(1:length(x))
  ##Combine two vectors as pairs
  mypair <- cbind(x,y,poi)
  ##get C Matrix
  Cmatrix <- cmat(x,y)
  ##Set initial values to rearrangement 
  Rearrange <- NULL
  ##keep recording the number of columns of C Matrix
  n <- ncol(mypair)*length(Cmatrix)/length(mypair)
  
  while(n >1){
    ##get row sums
    rowsums <- apply(X=Cmatrix,MARGIN=1,FUN=sum)
    ## A break situation
    if(max(rowsums)-min(rowsums)==0) {break}
    ##record the index of the row with smallest sum
    MinIndex <- order(rowsums)[1]
    ##Rearrange the pairs
    Rearrange <- rbind(Rearrange,mypair[MinIndex,])  
    ##Renew Cmatrix
    Cmatrix <- Cmatrix[-MinIndex,-MinIndex]
    ##Renew mypair
    mypair <- mypair[-MinIndex,]
    n <- ncol(mypair)*length(Cmatrix)/length(mypair)
  }
  ##Reverse the output
  Rearrange <- cbind(rev(Rearrange[,1]),rev(Rearrange[,2]),rev(Rearrange[,3]))
  Rearrange <- rbind(mypair,Rearrange)
  rownames(Rearrange) <- NULL
  return(Rearrange)
}

x<-Abnd_Measure
y<-Prev_Measure

Reorder.transp <- ReOrder(x,y)
x2 <- Reorder.transp[,1]
y2 <- Reorder.transp[,2]
t2 <- TauPath(cmat(x2,y2))
#plot(t2)

k <- 50
sig <- Reorder.transp[1:k,1:3]
index <- sig[,3]
rank <- c(1:length(Abnd_Measure))
strng_index<- rank[index]
weak_index <- rank[-index]
a<-strng_index

#Selecting the top K features from tau-path method
Data1<-Data[,a_int]

#Aggregating the remaining features into a single feature
Data1$Others<-rowSums(Data[,-a])
pred_name=colnames(Data1)


topic1<-c("Topic1","Topic2","Topic3")



#class labels associated with the data
load(file = "actual.rda")
summary(actual)
BarChart(actual,ylab =" Count-Health Status",xlab="Health Status")


data<-cbind(Data,actual)
#Boruta Method

set.seed(3)
boruta.bank_train <- Boruta(actual~., data = data, doTrace = 2,maxruns=1000)
print(boruta.bank_train)
col1<-c("Clostridium.sensu.stricto", "Eggerthella", "Paraprevotella",
        "unclassified_Proteobacteria")
data=Data1[,col1]
data$Others<-c()
c11<-which(colnames(Data) %in% col1)
data$Others<-rowSums(Data1[,-c11])
#rowSums(data)
library(tm)
head(data)
myDTM<- as.DocumentTermMatrix(data, weighting = weightTf)
myDTM
library(topicmodels)
lda_mod<- LDA(myDTM, k = 3,control = list(seed=3))
Gamma_lda<-lda_mod@gamma
pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
#for k topics k! combinations.
library(combinat)
classcomb<-permn(3)
accuracy1<-c()
WAcc1<-function(tab)
{
  a1<-((0.6*(diag(tab)[1]/table(actual)[1]))+(0.15*(diag(tab)[2]/table(actual)[2]))+(0.25*(diag(tab)[2]/table(actual)[2])))
  return(as.numeric(a1))
}
for(i in 1:length(classcomb))
{
  pred11<-factor(pred,levels=topic1[classcomb[[i]]])
  tab1<-table(pred11,actual);
  accuracy1[i]<-WAcc1(tab1)
}
b1=which.max(accuracy1)
pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
tab1<-table(pred11,actual);tab1
#sum(diag(tab1))/sum(tab1)
qf<-WAcc1(tab1)
round(qf*100,2)

