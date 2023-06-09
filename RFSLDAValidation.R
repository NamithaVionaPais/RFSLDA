setwd("/Users/namithapais/Documents/Documents - Namitha’s MacBook Pro/Thesis_NVP/Chapter1/Paper_Github/")

load(file = "MicrobiomeData_109.rda")
head(Data_109)
col1<-c("Response","Clostridium.XlVa","unclassified_Bacteroidales","Paraprevotella","Dialister","Coprococcus","Enterobacter",
        "unclassified_Fusobacteriaceae","Gemella","unclassified_Proteobacteria","Comamonas" )         
#Split into train test
p <- 0.8
strats <- Data_109$Response
rr <- split(1:length(strats), strats)
set.seed(3)
idx <- sort(as.numeric(unlist(sapply(rr, function(x) sample(x, length(x) * p)))))

train <- Data_109[idx, ]

test <- Data_109[-idx, ]

summary(Data_109$Response)/nrow(Data_109)
summary(train$Response)/nrow(train)
summary(test$Response)/nrow(test)
actual<-test$Response
WAcc1<-function(tab)
{
  a1<-((0.6*(diag(tab)[1]/table(actual)[1]))+(0.15*(diag(tab)[2]/table(actual)[2]))+(0.25*(diag(tab)[2]/table(actual)[2])))
  return(as.numeric(a1))
}

#Multinomial Logit

library(nnet)

# Training the multinomial model
multinom_model <- multinom(Response ~ ., data = train)
summary(multinom_model)


# Predicting the values for train dataset
prediction_mlogit_train <- predict(multinom_model, newdata = train, "class")
# Building classification table
mlogit_tab_train <- table(prediction_mlogit_train,train$Response);mlogit_tab_train
# Calculating accuracy - sum of diagonal elements divided by total obs
round((sum(diag(mlogit_tab_train))/sum(mlogit_tab_train)),2)

# Predicting the class for test dataset
prediction_mlogit_test <- predict(multinom_model, newdata = test, "class")
# Building classification table
mlogit_tab_test<- table(prediction_mlogit_test,test$Response);mlogit_tab_test
mlogit_tab_test
round((sum(diag(mlogit_tab_test))/sum(mlogit_tab_test)),2)
WAcc1(mlogit_tab_test)
#MMulticlass SVM
library(e1071) 

svm_model<- svm(Response~., data=train, 
                method="C-classification", kernal="linear", 
                gamma=0.1, cost=10)

summary(svm_model)

prediction_svm_train <- predict(svm_model, train)
svm_tab_train<- table(train$Response, prediction_svm_train);svm_tab_train

round((sum(diag(svm_tab_train))/sum(svm_tab_train)),2)
prediction_svm_test <- predict(svm_model, test)
svm_tab_test<- table(prediction_svm_test,test$Response);svm_tab_test
svm_tab_test
round((sum(diag(svm_tab_test))/sum(svm_tab_test)),2)
WAcc1(svm_tab_test)

library(topicmodels)
train <- Data_109[idx, col1]
test <- Data_109[-idx, col1]
train1<-train[,-1]
myDTM<- as.DocumentTermMatrix(train1, weighting = weightTf)
myDTM
lda_mod<- LDA(myDTM, k = 3,control = list(seed=33))
Gamma_lda<-lda_mod@gamma
topic1<-c("Topic1","Topic2","Topic3")
pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
#for k topics k! combinations.
classcomb<-permn(3)
#list
accuracy1<-c()
actual<-train[,1]
for(i in 1:length(classcomb))
{
  pred11<-factor(pred,levels=topic1[classcomb[[i]]])
  tab1<-table(pred11,train[,1]);
  accuracy1[i]<-WAcc1(tab1)
}
b1=which.max(accuracy1)
pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
tab1<-table(pred11,train[,1]);tab1
sum(diag(tab1))/sum(tab1)
lda_inf <- posterior(lda_mod, test[,-1])
pred_test<-as.factor(topic1[(apply(lda_inf$topics, 1, which.max))])
t1<-c()
for(i in 1:3)
{
  t1[i]<-paste("Topic",unlist(classcomb[b1])[i],sep="")
}

pred1_test1<-c()
for( i in 1:length(pred_test))
{
  if(pred_test[i]==t1[1]){
    pred1_test1[i]="Healthy"
  }else if(pred_test[i]==t1[2]){
    pred1_test1[i]="Infection"
  }else if(pred_test[i]==t1[3]){
    pred1_test1[i]="Stress"
  }
}

actual<-test$Response
tabf<-table(pred1_test1,test[,1]);tabf
sum(diag(tabf))/length(test[,1])
WAcc1(tabf)

