Data_sub<-Data
tot_sub<-rowSums(Data_sub)
#Calculating Abundance
Abnd_Data_sub<-Data_sub
for(i in 1:nrow(Data_sub))
{
  Abnd_Data_sub[i,]<-(Data_sub[i,]/tot_sub[i])*100
}

Abnd_Measure<-colSums(Abnd_Data_sub)
th<-1

#Calculating Prevelance
Prev_Data_sub<-ifelse(Abnd_Data_sub<=1,0,1)
Prev_Measure<-colSums(Prev_Data_sub)


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
union(a1,union(a2,union(a3,a)))
length(union(a1,union(a2,union(a3,a))))

a_int<-union(a1,union(a2,union(a3,a)))
