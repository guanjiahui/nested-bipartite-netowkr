library(bipartite)
load("ecology.RData")
#checker-board
Swap=function(Adj,i,j,i1,j1){
  swap=0
  if(Adj[i,j]==0 & Adj[i1,j]==1 & Adj[i,j1]==1 & Adj[i1,j1]==0){
    Adj[i,j]=1 
    Adj[i1,j]=0 
    Adj[i,j1]=0 
    Adj[i1,j1]=1
    swap=1
    return(list(Adj,swap))
  }
  if(Adj[i,j]==1 & Adj[i1,j]==0 & Adj[i,j1]==0 & Adj[i1,j1]==1){
    Adj[i,j]=0 
    Adj[i1,j]=1 
    Adj[i,j1]=1 
    Adj[i1,j1]=0
    swap=1
    return(list(Adj,swap))
  }# end if
  else
    return(list(Adj=Adj,swap=swap))
}#Swap
############
CheckerBoard=function(Adj,io,jo){
  for (i in 1: (nrow(Adj)-1)){
    for (j in 1:(ncol(Adj)-1)){
      if(i!=io & j!=jo){
        if (Swap(Adj,i,j)[[2]]==1){
          Matrix=Swap(Adj,i,j)[[1]]
          i1=i
          j1=j
          break
        }
      }
    }#end for j
  }#end for i
  # Energy=GetBipEnergy(Matrix)
  #Temperature=nestedtemp(Matrix)$statistic
  # NODF=nestednodf(Matrix,order=FALSE)[3]$statistic[3]
  return(list(Matrix=Matrix,i1=i1,j1=j1))#,Energy=Energy,Temperature=Temperature, NODF=NODF))
}
###########'
###
CheckerBoard=function(Adj){
  i=1;j=1;i1=2;j1=2
  while(Swap(Adj,i,j,i1,j1)[[2]]!=1){
    A=sample(nrow(Adj),2,replace=TRUE)
    B=sample(ncol(Adj),2,replace=TRUE)
    i=A[1]
    j=B[1]
    i1=A[2]
    j1=B[2]
    Matrix=Swap(Adj,i,j,i1,j1)[[1]]
  }
  # Energy=GetBipEnergy(Matrix)
  #Temperature=nestedtemp(Matrix)$statistic
  # NODF=nestednodf(Matrix,order=FALSE)[3]$statistic[3]
  return(Matrix=Matrix)#,Energy=Energy,Temperature=Temperature, NODF=NODF))
}
###########
Checker=function(Adj,swap=1,iter,Mat=NULL){
  if(is.null(Mat))
    Mat<- list()
  
  if (swap==iter)
    return(Mat)
  #Stat=CheckerBoard(Adj,io,jo)
  mat=CheckerBoard(Adj)
  Mat[[swap]]=mat
  #i1=Stat$i1
  #j1=Stat$j1
  Mat<- Checker(mat,swap+1,iter,Mat)
  return(Mat)
}#Checker
####
###########
CheckerEnergy=function(Adj,swap=1,iter,Energy=NULL){
  if(is.null(Energy))
    Energy<- c()
  
  if (swap==iter)
    return(Energy)
  #Stat=CheckerBoard(Adj,io,jo)
  mat=CheckerBoard(Adj)
  Energy[swap]=GetBipEnergy(mat)
  #i1=Stat$i1
  #j1=Stat$j1
  Energy<- CheckerEnergy(mat,swap+1,iter,Energy)
  return(Energy)
}#Checker
iter=1000
EnergyChange=CheckerEnergy(mam3,swap=1,iter)
#######
#try using for loop 
iter=1000
temp=list()
temp[[1]]=CheckerBoard(mam3)
for (i in 2:iter){
  temp[[i]]=CheckerBoard(temp[[(i-1)]])
}
  
Etemp=sapply(1:iter,function(i) GetBipEnergy(temp[[i]]))





#######################


iter=1000
BinaryChange=Checker(mam2,swap=1,iter)
Echange=sapply(1:(iter-1),function(i) GetBipEnergy(BinaryChange[[i]]))
NODFchange=sapply(1:(iter-1),function(i) nestednodf(BinaryChange[[i]],order=FALSE)[3]$statistic[3])
Tchange=sapply(1:(iter-1),function(i) nestedness(BinaryChange[[i]])$temperature)
save(BinaryChange,Echange,NODFchange,Tchange,file="checker.RData")
##########
# checkerboard for 2x2 blocks
iter=1000
BinaryM1=Checker(mam2[1:5,1:20],swap=1,iter)
BinaryM2=mam2[1:5,21:26]
BinaryM3=Checker(mam2[6:28,1:20],swap=1,iter)
BinaryM4=Checker(mam2[6:28,21:26],swap=1,iter)
BinaryM42=Checker(BinaryM4[[49]],swap=1,iter)


Binary22=rbind(cbind(BinaryM1,BinaryM2),cbind(BinaryM3,BinaryM4))
Echange=sapply(1:(iter-1),function(i) GetBipEnergy(Binary22[[i]]))
NODFchange=sapply(1:(iter-1),function(i) nestednodf(Binary22[[i]],order=FALSE)[3]$statistic[3])
Tchange=sapply(1:(iter-1),function(i) nestedness(Binary22[[i]])$temperature)
##########



#############
#function of N+ index(nestedness index from Patterson (1986))
mam4=t(mam3)[rev(1:ncol(mam3)),]
mamOrder=data.matrix(mamm[,2:ncol(mamm)])
########
Nplus <- function(Adj){
  m=nrow(Adj)
  n=ncol(Adj)
  q=colSums(Adj)
  c=numeric(m)
  for (i in 1:m){
    for (j in 1:n){
      if (Adj[i,j]==1)
        c[i]=min(q[j])
    }
  }
  temp=0
  for (i in 1:m){
    for (j in 1:n){
      if (Adj[i,j]==0 & q[j]>c[i])
        temp=temp+1
    }
  }#end for loop
  return(temp)
}#function Nplus
#######
Temp=lapply(1:50000, function(i) t(temp[[i]])[rev(1:ncol(temp[[i]])),])
Nchange=sapply(1:50000,function(i) Nplus(Temp[[i]]))
######
#check the plots 
plot(Etemp,type="l",main="Etemp")
plot(NODFtemp,type="l",main="NODFtemp")
plot(Ttemp,type="l",main="Temperature")
plot(Nchange,type="l",main="N-Plus index")
plot(density(Nchange),main="N-plus index distribution")
######
iter=50000
NODFtemp2=sapply(1:iter,function(i) nestednodf(Temp[[i]],order=FALSE)[3]$statistic[3])
Ttemp2=sapply(1:iter,function(i) nestedtemp(Temp[[i]])$statistic)

#########
#subtract two matrix to see if checkerboard is connected 
Subtract1=Matrix22[[77]]-Matrix22[[22]]
Subtract2=Matrix22[[100]]-Matrix22[[167]]
Subtract3=Matrix22[[3]]-Matrix22[[123]]
Subtract4=Matrix22[[89]]-Matrix22[[26]]
Subtract5=Matrix22[[57]]-Matrix22[[98]]
Subtract6=Matrix22[[12]]-Matrix22[[180]]
#######
#function for the subtraction matrix
Swap2=function(Adj,i,j,i1,j1){
  swap=0
  if(Adj[i,j]==-1 & Adj[i1,j]==1 & Adj[i,j1]==1 & Adj[i1,j1]==-1){
    Adj[i,j]=0
    Adj[i1,j]=0 
    Adj[i,j1]=0
    Adj[i1,j1]=0
    swap=1
    return(list(Adj,swap))
  }
  if(Adj[i,j]==1 & Adj[i1,j]==-1 & Adj[i,j1]==-1 & Adj[i1,j1]==1){
    Adj[i,j]=0
    Adj[i1,j]=0 
    Adj[i,j1]=0 
    Adj[i1,j1]=0
    swap=1
    return(list(Adj,swap))
  }# end if
  else
    return(list(Adj=Adj,swap=swap))
}#Swap2
######
CheckerBoard2=function(Adj){
  i=1;j=1;i1=2;j1=2
  while(Swap2(Adj,i,j,i1,j1)[[2]]!=1){
    A=sample(nrow(Adj),2,replace=TRUE)
    B=sample(ncol(Adj),2,replace=TRUE)
    i=A[1]
    j=B[1]
    i1=A[2]
    j1=B[2]
    Matrix=Swap2(Adj,i,j,i1,j1)[[1]]
  }
  # Energy=GetBipEnergy(Matrix)
  #Temperature=nestedtemp(Matrix)$statistic
  # NODF=nestednodf(Matrix,order=FALSE)[3]$statistic[3]
  return(Matrix=Matrix)#,Energy=Energy,Temperature=Temperature, NODF=NODF))
}
########
Checker2=function(Adj,swap=1,Mat=NULL){

  if(is.null(Mat))
    Mat<- list()


  mat=CheckerBoard2(Adj)
  Mat[[swap]]=mat

  if (sum(Mat[[swap]]==0)==length(Mat))
    return(Mat)
  
  Mat<- Checker2(mat,swap+1,Mat)
  return(Mat)
  

}#Checker
########
CheckerBoard3=function(Adj){
  indi=which(Adj==1,arr.ind=T)
  indj=which(Adj==-1,arr.ind=T)
  i=1;j=1;i1=2;j1=2
  while(Swap2(Adj,i,j,i1,j1)[[2]]!=1){
    A=sample(nrow(indi),2,replace=TRUE)
    i=indi[A[1],][1]
    j=indi[A[1],][2]
    i1=indi[A[2],][1]
    j1=indi[A[2],][2]
  }#while loop   
  Matrix=Swap2(Adj,i,j,i1,j1)[[1]]
  return(Matrix=Matrix)
}#CheckerBoard3()
#######
Sub1=Checker2(Subtract1,swap=1)
a=sum(Subtract1==-1)
iter=10
Sub1=list()
Sub1[[1]]=Subtract1
for (i in 2:iter){
  Sub1[[i]]=CheckerBoard2(Sub1[[(i-1)]])
}
for (i in 13:16)
  Sub1[[i]]=CheckerBoard3(Sub1[[(i-1)]])
####
Sub2=list()
Sub2[[1]]=Subtract2
for (i in 2:18){
  Sub2[[i]]=CheckerBoard3(Sub2[[(i-1)]])
  print(i)
}
##################
#recurison time
rec=c()
for (i in 2:length(ind1)){
  temp=ind1[i]-ind1[i-1]
  rec[i]=temp
}
ind.new=which(rec!=1)
recur=rec[ind.new]
####
rec2=c()
ind2=which(Nchange<67)
for (i in 2:length(ind2)){
  temp=ind2[i]-ind2[i-1]
  rec2[i]=temp
}
ind.new2=which(rec2!=1)
recur2=rec2[ind.new2]
##
rec3=c()
ind3=which(Ttemp>quantile(Ttemp,probs=0.9))
for (i in 2:length(ind3)){
  temp=ind3[i]-ind3[i-1]
  rec3[i]=temp
}
ind.new3=which(rec3!=1)
recur3=rec3[ind.new3]
###################