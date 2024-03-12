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
Data1$Others<-rowSums(Data[,-a_int])
pred_name=colnames(Data1)


topic1<-c("Topic1","Topic2","Topic3")



########################## Reference Model ########################
#Model when all the bacteria types are used as features 
data=Data1
myDTM<- as.DocumentTermMatrix(data, weighting = weightTf)
myDTM
lda_mod<- LDA(myDTM, k = 3,control = list(seed=6))
Gamma_lda<-lda_mod@gamma
pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
WAcc1<-function(tab)
{
  a1<-((0.6*(diag(tab)[1]/table(actual)[1]))+(0.15*(diag(tab)[2]/table(actual)[2]))+(0.25*(diag(tab)[3]/table(actual)[3])))
  return(as.numeric(a1))
}
#for k topics k! combinations.
classcomb<-permn(3)
accuracy1<-c()
for(i in 1:length(classcomb))
{
  pred11<-factor(pred,levels=topic1[classcomb[[i]]])
  tab1<-table(pred11,actual);
  accuracy1[i]<-WAcc1(tab1)
}
b1=which.max(accuracy1);b1
pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
tab1<-table(pred11,actual);tab1

qf<-WAcc1(tab1);qf
round(qf,4)
rownames(Gamma_lda)<-rownames(Data1)


########################## RFSLDA ########################
#hyperparamters
niter1=1
q1=1
q2=0
inf_prop=0
str_prop=0
health_prop=0
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
  lda_mod<- LDA(myDTM, k = 3,control = list(seed=3))
  Gamma_lda<-lda_mod@gamma
  pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
  #for k topics k! combinations.
  classcomb<-permn(3)
  #list
  accuracy1<-c()
  
  for(i in 1:length(classcomb))
  {
    pred11<-factor(pred,levels=topic1[classcomb[[i]]])
    tab1<-table(pred11,actual);
    accuracy1[i]<-WAcc1(tab1)
  }
  b1=which.max(accuracy1)
  pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
  tab1<-table(pred11,actual);tab1
  qint<-WAcc1(tab1);qint
  qf1<-qint
  #while(abs(q1-q2)>threshold|health_prop<0.7| inf_prop<0.5|str_prop<0.5)
  while(abs(qf1-qf)>threshold|qf<0.7)
  {
    qf1<-qf
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
    lda_mod1<- LDA(myDTM1, k = 3,control = list(seed=3))
    Gamma_lda1<-lda_mod1@gamma
    pred1<-as.factor(topic1[(apply(Gamma_lda1, 1, which.max))])
    #for k topics k! combinations.
    classcomb<-permn(3)
    #list
    accuracy1<-c()
    for(i in 1:length(classcomb))
    {
      pred11<-factor(pred1,levels=topic1[classcomb[[i]]])
      tab1<-table(pred11,actual);
      accuracy1[i]<-WAcc1(tab1)
    }
    b1=which.max(accuracy1)
    pred11<-factor(pred1,levels=topic1[classcomb[[b1]]])
    tab1<-table(pred11,actual);tab1
    q1<-WAcc1(tab1)
    health_prop1=tab1[1,1]/table(actual)[1]
    inf_prop1=tab1[2,2]/table(actual)[2]
    str_prop1=tab1[3,3]/table(actual)[3]
    #Model 2
    myDTM2 <- as.DocumentTermMatrix(data2, weighting = weightTf)
    myDTM2
    lda_mod2<-LDA(myDTM2, k = 3,control = list(seed=3))
    Gamma_lda2<-lda_mod2@gamma
    pred2<-as.factor(topic1[(apply(Gamma_lda2, 1, which.max))])
    classcomb<-permn(3)
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
    q2<-WAcc1(tab2)
    health_prop2=tab2[1,1]/table(actual)[1]
    inf_prop2=tab2[2,2]/table(actual)[2]
    str_prop2=tab2[3,3]/table(actual)[3]
    if(q2>q1)
    {
      s1=s2
      health_prop=health_prop2
      inf_prop=inf_prop2
      str_prop=str_prop2
      qf=q2
    }else if(q2<=q1){
      u=exp(-c*(q1-q2))
      rb1=rbinom(1, 1, u)
      if(rb1==1)
      {
        s1=s1
        qf=q1
        health_prop=health_prop1
        inf_prop=inf_prop1
        str_prop=str_prop1
      }else if(rb1==0){
        s1=s2
        qf=q2
        health_prop=health_prop2
        inf_prop=inf_prop2
        str_prop=str_prop2
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
  lda_mod<- LDA(myDTM, k = 3,control = list(seed=3))
  Gamma_lda<-lda_mod@gamma
  pred<-as.factor(topic1[(apply(Gamma_lda, 1, which.max))])
  #for k topics k! combinations.
  classcomb<-permn(3)
  #list
  accuracy1<-c()
  for(i in 1:length(classcomb))
  {
    pred11<-factor(pred,levels=topic1[classcomb[[i]]])
    tab1<-table(pred11,actual);
    accuracy1[i]<-WAcc1(tab1)
  }
  b1=which.max(accuracy1)
  pred11<-factor(pred,levels=topic1[classcomb[[b1]]])
  tab1<-table(pred11,actual);tab1
  qf<-WAcc1(tab1)
  print("Done")
  return(list(qint,pred_set1,qf,niter1,colnames(data),tab1,sf))
}


########################## Optimal Set ########################
#n=number of features available in the predictor basket
n=length(colnames(Data1))
#p=number of predictors to initialize the random feature set
p=(n/4)
#niter-Minimum number of iterations to be performed
niter=5000
#threshold=threshold for difference in accuracy
threshold=0.000001
c=1

get_lda_results<-function(seed)
{
  set.seed(seed)
  niter1=1
  q1=1
  q2=0
  health_prop=0
  inf_prop=0
  str_prop=0
  qf=0
  return(lda_fs(n,p,niter,threshold,c))
}

lda_final_iters<-list()
R<-50
for(i in 1:R)
{
  lda_final_iters[[i]]<-get_lda_results(i)
}

#save(lda_final_iters, file="lda_final_iters_fromtop80.RData") 

#From the R iterations identify the iteration with the highest weighted accuracy
R<-50
#load("lda_final_iters_fromtop50.RData") 
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
a12=unlist(lda_final_iters[[(which.max(unlist(ac_fin_list)))]][[5]])
which.max(ac_fin_list)
data=Data1[,a12]

#Final Model
library(tm)
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
  a1<-((0.6*(diag(tab)[1]/table(actual)[1]))+(0.15*(diag(tab)[2]/table(actual)[2]))+(0.25*(diag(tab)[3]/table(actual)[3])))
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

########################## Visualization  ########################
Data_result<-data.frame(matrix(nrow=89,ncol=6))
colnames(Data_result)<-c("Subject_ID","Infection_Prop","Healthy_Prop","Stress_Prop","Predicted_Cluster","CL4")
Data_result$Subject_ID<-row.names(Data)
Data_result[,2:4]<-round(Gamma_lda,4)
Data_result$Predicted_Cluster<-pred11
Data_result$CL4<-actual
#View(Data_result)
#write.csv(Data_result,"Data_result_subject_specific.csv")

#Visualizing clusters via composition proportions
colnames(Gamma_lda)=c("Cluster 1","Cluster 2","Cluster 3")
#colnames(Gamma_lda)=c("Stress","Healthy","Infection")
library(DirichletReg)
plot(DR_data(Gamma_lda), a2d = list(colored = TRUE, c.grid = TRUE, col.scheme = c("dims")))

#Top bacteria across different health status with their estimated proportions
library(tidytext)
ap_topics <- tidy(lda_mod, matrix = "beta")
ap_topics$HS<-NA
ap_topics$HS<-ifelse(ap_topics$topic==1,"Stress",ap_topics$HS)
ap_topics$HS<-ifelse(ap_topics$topic==2,"Infection",ap_topics$HS)
ap_topics$HS<-ifelse(ap_topics$topic==3,"Healthy",ap_topics$HS)
library(ggplot2)

library(dplyr)
ap_top_terms <- ap_topics %>%
  group_by(topic) %>%
  slice_max(beta, n = 10) %>% 
  ungroup() %>%
  arrange(topic, -beta)
ap_top_terms %>%
  mutate(term = reorder_within(term, beta, topic)) %>%
  ggplot(aes(beta, term, fill = factor(topic))) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ HS, scales = "free_y") +
  scale_y_reordered()

#Composition of optimal bacteria types across subjects based on their read counts

library(ggplot2)
library(viridis)
library(hrbrthemes)
library(reshape2)
data_long <- data/rowSums(data)             
data_long$subgroup <- as.factor(rownames(data_long))
data_long <- melt(data_long, id.vars = "subgroup")
data_long 
data_long$Bacteria_Types=data_long$variable

ggplot(data_long,                  # Stacked barplot using ggplot2
       aes(x = subgroup,
           y = value,
           fill = Bacteria_Types)) +
  geom_bar(stat = "identity")+ labs(y = "Read Counts")+
  labs(x = "Subjects")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank() 
  )

#Within-group similarity across each health status level using the t-SNE method
library(Rtsne)
library(ggplot2)
library(ggpubr)
head(Gamma_lda)
cl=c()
for(i in 1:nrow(Gamma_lda))
{
  if(max(Gamma_lda[i,]>0.7))
  {
    cl[i]<-which.max(Gamma_lda[i,])
  } else {
    cl[i]<-paste(which(rank(-Gamma_lda[i,])==1),which(rank(-Gamma_lda[i,])==2),sep="-")
  }
  
}
cl<-as.factor(cl)

cl1<-c()
cl1[which(cl==1)]="Predominantly-Stress"
cl1[which(cl==3)]="Predominantly-Healthy"
cl1[which(cl==2)]="Predominantly-Infection"
cl1[which(cl=="1-3")]="Stress-Healthy"
cl1[which(cl=="1-2")]="Stress-Infection"
cl1[which(cl=="3-1")]="Healthy-Stress"
cl1[which(cl=="3-2")]="Healthy-Infection"
cl1[which(cl=="2-1")]="Infection-Stress"
cl1[which(cl=="2-3")]="Infection-Healthy"

cl1<-as.factor(cl1)
cl1 <- ordered(cl1, levels = c("Predominantly-Healthy","Healthy-Infection","Healthy-Stress",
                               "Predominantly-Infection","Infection-Healthy","Infection-Stress",
                               "Predominantly-Stress","Stress-Healthy","Stress-Infection"))
summary(cl1)

pred1<-c()
pred1[which(pred=="Topic1")]="Stress"
pred1[which(pred=="Topic2")]="Infection"
pred1[which(pred=="Topic3")]="Healthy"
pred1<-as.factor(pred1)

# Calculating T-SNE's Y values which indicate where each cluster is on the 2-D map
set.seed(3)
y_tsne <- Rtsne(Gamma_lda,pca=FALSE,perplexity=29, max_iter=1000, check_duplicates = FALSE) # Run TSNE

# Create a data frame using IRIS values and the T-SNE values
df_to_gg<-as.data.frame(cbind(cl1, as.data.frame(y_tsne$Y)))

Health_Groups<-pred1
# Specify column names
names(df_to_gg)<-c("Similarity", "Y.1", "Y.2")
library(RColorBrewer)
# Show the objects in the 2D T-SNE representation
ggplot(df_to_gg, aes(x=Y.1, y=Y.2, color=Similarity,pch=(Health_Groups)))+geom_point()+theme_minimal() +
  labs(title=paste("" )) +
  scale_color_manual(values = brewer.pal(n = 9, name = "Set1"))+guides(color = guide_legend(
    override.aes=list(shape = 18)))





