source("Bootstrapping (1).r")
bootstrap.energy=function(rowcut,colcut,rowTree,colTree,Binary,iter){
  m=rowcut
  n=colcut
  group1=cutree(rowTree,k=m)
  group2=cutree(colTree,k=n)
  adj=list()
  k=1
  for(i in 1:m){
    for (j in 1:n){
      adj[[k]]=Binary[which(group1==i),which(group2==j)]
      k=k+1
    }
  }#end of for loop
  
  Energy=numeric(iter)
  for (i in 1:iter){
    Energy1=sapply(1:(m*n), function(k) Bootbinary(adj[[k]])$Energy)
    Energy[i]=sum(Energy1)
  }
  return(Energy)
}#bootstrap.energy

load("DCG.RData")
nb=1000
overall.both=bootstrap.energy(1,1,phy.an.tree,phy.pl.tree,BEEH,nb)
overall.col=bootstrap.energy(1,1,phy.an.tree,DCG.1,BEEH,nb)
overall.row=bootstrap.energy(1,1,DCG.1.lx,phy.pl.tree,BEEH,nb)
E.bootc=numeric(nb)
BEEH.bootc=list()
for (i in 1:nb){
  a=Bootbinary(BEEH)
  E.bootc[i]=a$Energy
  BEEH.bootc[[i]]=a$Matrix
}
save(overall.both,overall.col,overall.row,E.bootc,BEEH.bootc,file="BEEH_energy.RData")
###############
par(mfrow=c(1,1))
plot(overall,ylim=c(-480,-210),main="BEEH",type="l")
points(BEEH.Ec,col="red")
#confidence interval(percentile)
CI=quantile(overall.c,probs=c(0.025,0.975))
##########
#hypothesis testing 
E.star=overall
E.null=E1.c
E.hat=GetBipEnergy(BEEH[DCG.10$order,DCG.11$order])
t.hat=quantile(abs(E.star-E.hat),prob=0.05)
test.stat=abs(E.hat-E.null)
print(c(t.hat,test.stat))
#####
#permutation
plot(E.permute,ylim=c(-480,max(E.permute)),main="BEEH",type="l")
points(BEEH.Ec,col="red")
CI.P=quantile(E.permute,probs=c(0.025,0.975))

#####