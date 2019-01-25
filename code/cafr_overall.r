source("bootstrapEnergy.r")
load("DCG.RData")
nb=800
cafrBootBoth=bootstrap.energy(1,1,phy.an.CAFR,phy.pl.CAFR,Binary_CAFR,nb)
cafrBootCol=bootstrap.energy(1,1,phy.an.CAFR,DCG_tree1,Binary_CAFR,nb)
cafrBootRow=bootstrap.energy(1,1,DCG_tree1.lx,phy.pl.CAFR,Binary_CAFR,nb)
cafrBootCoup=bootstrap.energy(1,1,DCG_tree10,DCG_tree9,Binary_CAFR,nb)
E.cafr=numeric(nb)
cafrBootMatrix=list()
for (i in 1:nb){
  a=Bootbinary(Binary_CAFR)
  E.cafr[i]=a$Energy
  cafrBootMatrix[[i]]=a$Matrix
}
save(cafrBootBoth,cafrBootCol,cafrBootRow,cafrBootCoup,E.cafr,cafrBootMatrix,file="CAFR_energy.RData")
dat.cafr=data.frame(dens=c(cafrBootCol,cafrBootRow,cafrBootCoup,E.cafr),lines=rep(c("col","row","coup","initial"),each=nb))
densityplot(~dens,data=dat.cafr,groups = lines,plot.points = FALSE, ref = TRUE,
            auto.key=list(lines=TRUE),main="CAFR bootstrap both col phylogeny")
##################################
E1.cafr=GetBipEnergy(Binary_CAFR[DCG_tree10$order,DCG_tree9$order])
E2.cafr=GetBipEnergy(Binary_CAFR[phy.an.CAFR$order,DCG_tree1$order])
E3.cafr=GetBipEnergy(Binary_CAFR[DCG_tree1.lx$order,phy.pl.CAFR$order])
E4.cafr=GetBipEnergy(Binary_CAFR[phy.an.CAFR$order,phy.pl.CAFR$order])
Efinal.cafr=GetBipEnergy(Binary_CAFR[DCG_tree12$order,DCG_tree11$order])
cafrPoint=c(E1.cafr,E2.cafr,E3.cafr,E4.cafr)

plot(E.cafr,ylim=c(E1.cafr,max(E.cafr)),main="CAFR")
points(cafrPoint,col="green")
###################################
#hypothesis testing 
E.star=E.cafr
E.null=E4.cafr
E.hat=GetBipEnergy(Binary_CAFR)
t.hat=quantile(abs(E.star-E.hat),prob=0.05)
test.stat=abs(E.hat-E.null)
print(c(t.hat,test.stat))
########
#confidence interval 
CI_cafr=quantile(cafrBootRow,probs=c(0.025,0.975))
############
#permutation
n.random=5000
CAFR.permute=array(0,c(nrow(Binary_CAFR),ncol(Binary_CAFR),n.random))
for (i in 1:n.random){
  row.permute=sample(nrow(Binary_CAFR),replace=FALSE)
  col.permute=sample(ncol(Binary_CAFR),replace=FALSE)
  CAFR.permute[,,i]=Binary_CAFR[row.permute,]
  CAFR.permute[,,i]=CAFR.permute[,col.permute,i]
}
E.per.cafr=sapply(1:n.random, function(i) GetBipEnergy(CAFR.permute[,,i]))
CI_cafr.per=quantile(E.per.cafr,probs=c(0.025,0.975))
plot(E.per.cafr,ylim=c(min(cafrPoint),max(E.per.cafr)),main="CAFR permutation",type="l")
points(cafrPoint,col="red")