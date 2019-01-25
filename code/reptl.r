source("Bootstrapping (1).r")
source("EnergyOptim.r")
load("DCG.RData")
load("CoupGeo.RData")
###
iter=500
Energy.coarse=numeric(iter)
Energy.fine=numeric(iter)
for (l in 1:iter){
  sub1=Bootbinary(reptl[1:11,1:6])$Matrix
  sub2=Bootbinary(reptl[1:11,7:44])$Matrix  
  
  sub3=Bootbinary(reptl[12:20,1:6])$Matrix  
  sub4=Bootbinary(reptl[12:20,7:44])$Matrix
  
  Binarytotal=rbind(cbind(sub1,sub2),cbind(sub3,sub4))
  Binarytotal=cbind(Binarytotal,reptl[,45:47])
  Energy.coarse[l]=GetBipEnergy(Binarytotal)
}

for (l in 1:iter){
  sub11=Bootbinary(reptl[1:2,1:6])$Matrix
  sub12=Bootbinary(reptl[3:4,1:6])$Matrix
  sub13=Bootbinary(reptl[5:11,1:6])$Matrix
  sub14=Bootbinary(reptl[12:20,1:6])$Matrix
  
  sub21=Bootbinary(reptl[1:2,7:44])$Matrix
  sub22=Bootbinary(reptl[3:4,7:44])$Matrix
  sub23=Bootbinary(reptl[5:11,7:44])$Matrix
  sub24=Bootbinary(reptl[12:20,7:44])$Matrix
  Binarytotal=cbind(rbind(sub11,sub12,sub13,sub14),rbind(sub21,sub22,sub23,sub24))
  Binarytotal=cbind(Binarytotal,reptl[,45:47])
  Energy.fine[l]=GetBipEnergy(Binarytotal)
}
#########
save(Energy.coarse,Energy.fine,file="reptl.RData")
library(lattice)
dat.reptl<- data.frame(dens = c(Energy.coarse,Energy.fine,E.reptl.per), lines = rep(c("median","fine","coarse"), c(500,500,10000)))
densityplot(~dens,data=dat.reptl,groups = lines,plot.points = FALSE, ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="reptile")
###############################
iter=100
Energy.overall=numeric(iter)
for (l in 1:iter){
  Binary.no=Bootbinary(reptl)$Matrix
  Energy.noblock[l]=GetBipEnergy(Binary.no)
}
#######
plot(density(Energy.avif.fine))
####
n.random=10000
reptl.permute=array(0,c(nrow(reptl),ncol(reptl),n.random))
for (i in 1:n.random){
  row.permute=sample(nrow(reptl),replace=FALSE)
  col.permute=sample(ncol(reptl),replace=FALSE)
  reptl.permute[,,i]=reptl[row.permute,]
  reptl.permute[,,i]=reptl.permute[,col.permute,i]
}
E.reptl.per=sapply(1:n.random, function(i) GetBipEnergy(reptl.permute[,,i]))
########
n.random=10000
avif.permute=array(0,c(nrow(avif),ncol(avif),n.random))
for (i in 1:n.random){
  row.permute=sample(nrow(avif),replace=FALSE)
  col.permute=sample(ncol(avif),replace=FALSE)
  avif.permute[,,i]=avif[row.permute,]
  avif.permute[,,i]=avif.permute[,col.permute,i]
}
E.avif.per=sapply(1:n.random, function(i) GetBipEnergy(avif.permute[,,i]))
E.avif.per1=E.avif.per
########
dat.avif<- data.frame(dens = c(Energy.avif.fine,E.avif.per1), lines = rep(c("fine","coarse"), c(100,10000)))
densityplot(~dens,data=dat.avif,groups = lines,plot.points = FALSE, ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="Avifauna")