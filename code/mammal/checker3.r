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
  return(Matrix=Matrix)
}
###########

#try using for loop 
iter=1.0e+6
MatChecker=list()
MatChecker[[1]]=CheckerBoard(mam2)
H=c()
H[1]=159
for (i in 2:iter){
  MatChecker[[i]]=CheckerBoard(MatChecker[[(i-1)]])
  CC=array(MatChecker[[i]],dim=c(length(mam2)/4,4))
  H[i]=t(CC[,1])%*%CC[,2]%*%t(CC[,3])%*%CC[,4]
  if (i==250000)
    print("a quarter")
  if (i==5e+05)
    print("HALF")
}

save(MatChecker,H,file="checkerE10.RData")
#######
Hlist=unique(H)
indexH=lapply(1:length(Hlist), function(i) which(H==Hlist[i]))
length_index=sapply(1:length(Hlist), function(i) length(indexH[[i]]))

for (i in 1:length_index[1]){
  for(j in 1:length_index[1]){
    if (sum((MatChecker[[indexH[[1]][i]]]-MatChecker[[indexH[[1]][j]]])!=0)==0 &i!=j){
      Same1=i
      Same2=j
      print(paste(i,j))
      break
    }
  }
  if (sum((MatChecker[[indexH[[1]][i]]]-MatChecker[[indexH[[1]][j]]])!=0)==0 &i!=j)
    break
  print(i)
}

Same21=c(233397,233399)
Same11=c(389608,389610)
#########
cor.test(Etemp[300:50000],NODFtemp[300:50000],method="pearson",
         alternative=c("less"))
cor.test(Etemp[300:50000],Nchange[300:50000],method="pearson",
         alternative=c("less"))
cor.test(Etemp[300:50000],Ttemp[300:50000],method="pearson",
         alternative=c("two.sided"))
cor(Etemp,Nchange)
cor(Etemp,NODFtemp)
cor(Nchange,NODFtemp)
cor(Etemp,Ttemp)
cor(Ttemp,NODFtemp)
cor(Ttemp,Nchange)
#####
Matrix55=lapply(1:200,function(i) t(Tmam55[[i]])[rev(1:ncol(Tmam55[[i]])),])
Nplus55=sapply(1:200,function(i) Nplus(Matrix55[[i]]))
Energy55=sapply(1:200,function(i) GetBipEnergy(Tmam55[[i]]))
NODF55=sapply(1:200,function(i) nestednodf(Tmam55[[i]],order=FALSE)[3]$statistic[3])
Temperature55=sapply(1:200,function(i) nestedtemp(Tmam55[[i]])$statistic[3])
#####
EcheckerMost=sapply(1:1.0e+6, function(i) GetBipEnergy(MatChecker[[i]]))
