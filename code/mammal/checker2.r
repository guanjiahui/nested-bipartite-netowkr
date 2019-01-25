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

#try using for loop 
iter=50000
temp=list()
temp[[1]]=CheckerBoard(mam3)
for (i in 2:iter){
  temp[[i]]=CheckerBoard(temp[[(i-1)]])
}

Etemp=sapply(1:iter,function(i) GetBipEnergy(temp[[i]]))
NODFtemp=sapply(1:iter,function(i) nestednodf(temp[[i]],order=FALSE)[3]$statistic[3])
Ttemp=sapply(1:iter,function(i) nestedtemp(temp[[i]])$statistic)
save(temp,Etemp,NODFtemp,Ttemp,file="checker50.RData")