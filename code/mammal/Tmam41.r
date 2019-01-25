source("Bootstrapping (1).r")
load("ecology.RData")
library(bipartite)
##########
iter=200
Emam41=numeric(iter)
NODF41=list()
NODF41test=list()
Matrix41=list()
for (l in 1:iter){
  sub1=mam2[1:28,1:2]
  sub2=Bootbinary(mam2[1:28,3:11])$Matrix 
  sub3=Bootbinary(mam2[1:28,12:20])$Matrix
  sub4=Bootbinary(mam2[1:28,21:26])$Matrix

  Binarytotal=cbind(sub1,sub2,sub3,sub4)
  Tmam41[l]=nestedness(Binarytotal)$temperature
  NODF41[[l]]=nestednodf(Binarytotal,order=FALSE)
  NODF41test[[l]]=nestednodf(Binarytotal)
  Matrix41[[l]]=Binarytotal
}

save(Tmam41,NODF41,NODF41test,Matrix41,file="Tmammal41.RData")