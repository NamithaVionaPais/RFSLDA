setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/Thesis_NVP/Chapter1/NewData/")
load("cancer.RData")

dim(cancer)
head(cancer)
summary(cancer)

cancer$cls.sample<-as.factor(cancer$cls.sample)
p <- 0.8
strats <- cancer$cls.sample
rr <- split(1:length(strats), strats)
set.seed(3)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

train <- cancer[idx, ]

test <- cancer[-idx, ]

summary(train$cls.sample)/nrow(train)
summary(test$cls.sample)/nrow(test)
actual<-test$Response

#Multinomial Logit

library(nnet)

# Training the multinomial model
multinom_model <- multinom(cls.sample ~ ., data = train)
summary(multinom_model)


# Predicting the values for train dataset
prediction_mlogit_train <- predict(multinom_model, newdata = train, "class")
# Building classification table
mlogit_tab_train <- table(prediction_mlogit_train,train$cls.sample);mlogit_tab_train
# Calculating accuracy - sum of diagonal elements divided by total obs
round((sum(diag(mlogit_tab_train))/sum(mlogit_tab_train)),2)

# Predicting the class for test dataset
prediction_mlogit_test <- predict(multinom_model, newdata = test, "class")
# Building classification table
mlogit_tab_test<- table(prediction_mlogit_test,test$cls.sample);mlogit_tab_test
mlogit_tab_test
round((sum(diag(mlogit_tab_test))/sum(mlogit_tab_test)),2)

#MMulticlass SVM
library(e1071) 
svm_model<- svm(cls.sample~., data=train, 
                method="C-classification", kernal="linear", 
                gamma=0.1, cost=10)

summary(svm_model)

prediction_svm_train <- predict(svm_model, train)
svm_tab_train<- table(train$cls.sample, prediction_svm_train);svm_tab_train

round((sum(diag(svm_tab_train))/sum(svm_tab_train)),2)
prediction_svm_test <- predict(svm_model, test)
svm_tab_test<- table(prediction_svm_test,test$cls.sample);svm_tab_test
svm_tab_test
round((sum(diag(svm_tab_test))/sum(svm_tab_test)),2)

####RFSLDA
library(lessR)
library(tm)
library(topicmodels)
library(combinat)
myDTM<- as.DocumentTermMatrix(train[,-150], weighting = weightTf)
myDTM
lda_mod<- LDA(myDTM, k = 2,control = list(seed=33))
Gamma_lda<-lda_mod@gamma
topic1<-c("Topic1","Topic2")
library()
pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
#for k topics k! combinations.
classcomb<-permn(2)
#list
accuracy1<-c()
actual<-train[,150]
for(i in 1:length(classcomb))
{
  pred11<-factor(pred,levels=topic1[classcomb[[i]]])
  tab1<-table(pred11,train[,150]);
  accuracy1[i]<-sum(diag(tab1))/sum(tab1)
}
b1=which.max(accuracy1)
pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
tab1<-table(pred11,train[,150]);tab1
sum(diag(tab1))/sum(tab1)

#hyperparamters
Data1<-train[,-150]
pred_name=colnames(Data1)
niter1=1
q1=1
q2=0
qf=0
lda_fs<-function(n,p,niter,threshold,c)
{
  s1=sample(1:n,size=p,replace = FALSE)
  while(!all(rowSums(Data1[,s1])) > 0)
  {
    s1=sample(1:n,size=p,replace = FALSE)
  }
  
  pred_set1=pred_name[s1]
  sint=s1
  #Initial Model based on best set of features
  data=Data1[,sint]
  myDTM<- as.DocumentTermMatrix(data, weighting = weightTf)
  myDTM
  lda_mod<- LDA(myDTM, k = 2,control = list(seed=3))
  Gamma_lda<-lda_mod@gamma
  pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
  #for k topics k! combinations.
  classcomb<-permn(2)
  #list
  accuracy1<-c()
  
  for(i in 1:length(classcomb))
  {
    pred11<-factor(pred,levels=topic1[classcomb[[i]]])
    tab1<-table(pred11,actual);
    accuracy1[i]<-sum(diag(tab1))/sum(tab1)
  }
  b1=which.max(accuracy1)
  pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
  tab1<-table(pred11,actual);tab1
  qint<-sum(diag(tab1))/sum(tab1);qint
  while(abs(q1-q2)>threshold|qf<0.65)
  {
    #Tossing a coin
    #Total number of predictors in the basket
    n=length(colnames(Data1))
    #number of predictors to select the first pont
    coin=sample(c(1:3), size=1)
    if(coin==1)
    {
      a1=sample(s1,size=1)
      a2=sample((1:n)[-s1],size=1)
      s2=c(setdiff(s1,a1),a2)
    } else if(coin==2) {
      a1=sample(s1,size=1)
      s2=c(setdiff(s1,a1))
    } else if(coin==3) {
      a2=sample((1:n)[-s1],size=1)
      s2=c(s1,a2)
    }
    while(!all(rowSums(Data1[,s2])) > 0|length(s2)<3|length(s2)>(n-1))
    {
      coin=sample(c(1:3), size=1)
      if(coin==1)
      {
        a1=sample(s1,size=1)
        a2=sample((1:n)[-s1],size=1)
        s2=c(setdiff(s1,a1),a2)
      } else if(coin==2) {
        a1=sample(s1,size=1)
        s2=c(setdiff(s1,a1))
      } else if(coin==3) {
        a2=sample((1:n)[-s1],size=1)
        s2=c(s1,a2)
      }
    }
    data1=Data1[,s1]
    data2=Data1[,s2]
    #Model 1
    myDTM1 <- as.DocumentTermMatrix(data1, weighting = weightTf)
    myDTM1
    lda_mod1<- LDA(myDTM1, k = 2,control = list(seed=3))
    Gamma_lda1<-lda_mod1@gamma
    pred1<-as.factor(topic1[(apply(Gamma_lda1, 1, which.max))])
    #for k topics k! combinations.
    classcomb<-permn(2)
    #list
    accuracy1<-c()
    for(i in 1:length(classcomb))
    {
      pred11<-factor(pred1,levels=topic1[classcomb[[i]]])
      tab1<-table(pred11,actual);
      accuracy1[i]<-sum(diag(tab1))/sum(tab1)
    }
    b1=which.max(accuracy1)
    pred11<-factor(pred1,levels=topic1[classcomb[[b1]]])
    tab1<-table(pred11,actual);tab1
    q1<-sum(diag(tab1))/sum(tab1)
    #Model 2
    myDTM2 <- as.DocumentTermMatrix(data2, weighting = weightTf)
    myDTM2
    lda_mod2<-LDA(myDTM2, k = 2,control = list(seed=3))
    Gamma_lda2<-lda_mod2@gamma
    pred2<-as.factor(topic1[(apply(Gamma_lda2, 1, which.max))])
    classcomb<-permn(2)
    #list
    accuracy2<-c()
    for(i in 1:length(classcomb))
    {
      pred21<-factor(pred2,levels=topic1[classcomb[[i]]])
      tab2<-table(pred21,actual);
      accuracy2[i]<-sum(diag(tab2))/length(actual)
    }
    b2=which.max(accuracy2)
    pred21<-factor(pred2,levels=topic1[classcomb[[b2]]])
    tab2<-table(pred21,actual);tab2
    q2<-sum(diag(tab2))/sum(tab2)
    if(q2>q1)
    {
      s1=s2
      qf=q2
    }else if(q2<=q1){
      u=exp(-c*(q1-q2))
      rb1=rbinom(1, 1, u)
      if(rb1==1)
      {
        s1=s1
        qf=q1
      }else if(rb1==0){
        s1=s2
        qf=q2
      }
    }
    niter1=niter1+1
    print(paste(niter1,round(qf,4),sep="-"))
  }
  sf=s1
  #Final Model based on best set of features
  data=Data1[,sf]
  #Reference Model
  myDTM<- as.DocumentTermMatrix(data, weighting = weightTf)
  myDTM
  lda_mod<- LDA(myDTM, k = 2,control = list(seed=3))
  Gamma_lda<-lda_mod@gamma
  pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
  #for k topics k! combinations.
  classcomb<-permn(2)
  #list
  accuracy1<-c()
  for(i in 1:length(classcomb))
  {
    pred11<-factor(pred,levels=topic1[classcomb[[i]]])
    tab1<-table(pred11,actual);
    accuracy1[i]<-sum(diag(tab1))/sum(tab1)
  }
  b1=which.max(accuracy1)
  pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
  tab1<-table(pred11,actual);tab1
  qf<-sum(diag(tab1))/sum(tab1)
  print("Done")
  return(list(qint,pred_set1,qf,niter1,colnames(data),tab1,sf))
}


########################## Optimal Set ########################
#n=number of features available in the predictor basket
n=length(colnames(train[,-150]))
#p=number of predictors to initialize the random feature set
p=(n/2)
#niter-Minimum number of iterations to be performed
niter=500
#threshold=threshold for difference in accuracy
threshold=0.000001
c=1

get_lda_results<-function(seed)
{
  set.seed(seed)
  niter1=1
  q1=1
  q2=0
  qf=0
  return(lda_fs(n,p,niter,threshold,c))
}

lda_final_iters<-list()
R<-50
for(i in 1:R)
{
  lda_final_iters[[i]]<-get_lda_results(i)
}

#save(lda_final_iters, file="lda_final_iters_cancerdata.RData") 

#From the R iterations identify the iteration with the highest weighted accuracy
R<-50
load("lda_final_iters_cancerdata.RData") 
ac_int_list<-c()
ac_fin_list<-c()
pred_int_list<-list()
pred_fin_list<-list()
niter_list<-c()
tab_list<-list()
len_finl<-c()
for(i in 1:R)
{
  ac_int_list[i]<-lda_final_iters[[i]][1]
  ac_fin_list[i]<-lda_final_iters[[i]][3]
  pred_int_list[[i]]<-lda_final_iters[[i]][2]
  pred_fin_list[[i]]<-lda_final_iters[[i]][5]
  niter_list[i]<-lda_final_iters[[i]][4]
  tab_list[[i]]<-lda_final_iters[[i]][6]
  len_finl[i]<-length(unlist(lda_final_iters[[i]][5]))
}
which.max(unlist(ac_fin_list))
#a12=unlist(lda_final_iters[[(which.max(unlist(ac_fin_list)))]][[5]])
a12<-unlist(lda_final_iters[[41]][[5]])
which.max(ac_fin_list)
data=Data1[,a12]


#posterior(lda_mod, test[,-1])
myDTM<- as.DocumentTermMatrix(data, weighting = weightTf)
myDTM
lda_mod<- LDA(myDTM, k = 2,control = list(seed=3))
Gamma_lda<-lda_mod@gamma
topic1<-c("Topic1","Topic2")
pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
#for k topics k! combinations.
classcomb<-permn(2)
#list
accuracy1<-c()
actual<-train[,150]
for(i in 1:length(classcomb))
{
  pred11<-factor(pred,levels=topic1[classcomb[[i]]])
  tab1<-table(pred11,train[,150]);
  accuracy1[i]<-sum(diag(tab1))/sum(tab1)
}
b1=which.max(accuracy1)
pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
tab1<-table(pred11,train[,150]);tab1
sum(diag(tab1))/sum(tab1)

Datatest<-test[,-150]
lda_inf <- posterior(lda_mod, Datatest[,a12])
pred_test<-as.factor(topic1[(apply(lda_inf$topics, 1, which.max))])
t1<-c()
for(i in 1:2)
{
  t1[i]<-paste("Topic",unlist(classcomb[b1])[i],sep="")
}

pred1_test1<-c()
for( i in 1:length(pred_test))
{
  if(pred_test[i]==t1[1]){
    pred1_test1[i]=1
  }else if(pred_test[i]==t1[2]){
    pred1_test1[i]=2
  }
}
pred1_test1<-as.factor(pred1_test1)
actual<-as.factor(test[,150])
tabf<-table(pred1_test1,actual);tabf
sum(diag(tabf))/length(test[,1])


