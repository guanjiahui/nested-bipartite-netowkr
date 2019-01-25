source("Bootstrapping (1).r")
bootstrap.energy=function(rowcut,colcut,rowTree,colTree,Binary,iter){
  m=rowcut
  n=colcut
  group1=cutree(rowTree,k=m)
  group2=cutree(colTree,k=n)
  row.order.cut=c()
  col.order.cut=c()
  row.order.cut[1]=0
  col.order.cut[1]=0
  a=b=0
  
  for (i in 2:(m+1)){
    temp=length(which(group1==(i-1)))
    a=a+temp
    row.order.cut[i]=a
  }
  for (j in 2:(n+1)){
    temp=length(which(group2==(j-1)))
    b=b+temp
    col.order.cut[j]=b
  }
  
  adj=list()
  k=1
  
  for(i in 1:m){
    for (j in 1:n){
      row=rowTree$order[(row.order.cut[i]+1):row.order.cut[(i+1)]]
      col=colTree$order[(col.order.cut[j]+1):col.order.cut[(j+1)]]
      adj[[k]]=Binary[row,col]
      k=k+1
    }
  }#end of for loop 
  
  Energy=numeric(iter)
  
  for (l in 1:iter){
    Binary_total=matrix(0,nrow(Binary),ncol(Binary))
    BinaryBoot=sapply(1:(m*n), function(k) Bootbinary(adj[[k]]))
    k=1
    
    for(i in 1:m){
      for (j in 1:n){
        row=rowTree$order[(row.order.cut[i]+1):row.order.cut[(i+1)]]
        col=colTree$order[(col.order.cut[j]+1):col.order.cut[(j+1)]]
        Binary_total[row,col]=BinaryBoot[,k]$Matrix
        k=k+1
      }
    }#end of for loop 
    
    Energy[l]=GetBipEnergy(Binary_total)
  }#end of for iter loop
  
  return(Energy)
}#bootstrap.energy

load("DCG.RData")
library(lattice)
n.b=500
b1c=bootstrap.energy(2,2,phy.an.tree,phy.pl.tree,BEEH,n.b)
b2c=bootstrap.energy(3,4,phy.an.tree,phy.pl.tree,BEEH,n.b)
b3c=bootstrap.energy(5,6,phy.an.tree,phy.pl.tree,BEEH,n.b)
c1c=bootstrap.energy(2,2,phy.an.tree,DCG.1,BEEH,n.b)
c2c=bootstrap.energy(3,4,phy.an.tree,DCG.1,BEEH,n.b)
c3c=bootstrap.energy(3,7,phy.an.tree,DCG.1,BEEH,n.b)
########
save(b1c,b2c,b3c,c1c,c2c,c3c,file="row_both2.Rdata")