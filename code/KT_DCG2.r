load("KT.RData")

row3.kt1=DCG.distance(KT3,DCG3.kt1,by.row=FALSE, k=5,method="euclidean",replicate=FALSE,r=5)

m1=nrow(KT)
n1=ncol(KT)
lx.kt3=matrix(0,n1,n1)
for (i in 1:m1){
  for(j in 1:m1)
    lx.kt3[i,j]=sum(KT3[i,]!=KT3[j,])
}
r3.kt1=row3.kt1+lx.kt3

temp= c(0.01,0.05,0.25,0.35,0.7,10)
Ens3.kt2_2=Eigen.plot(temp, selected.id <- c( 1,2,3,4,5,6),r3.kt1)
DCG3.kt2_2=DCGtree.plot(num.clusters.selected=c(1,2,3,5,5,5),"mam 1st row tree",Ens3.kt2_2,temp)

save(Ens3.kt2_2, DCG3.kt2_2,file="KT_DCG20162.RData")