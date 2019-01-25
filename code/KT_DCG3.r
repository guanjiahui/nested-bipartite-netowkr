load("KT_2017.RData")


m1=nrow(KT)
n1=ncol(KT)
lx.kt=matrix(0,m1,m1)
for (i in 1:m1){
  for(j in 1:m1)
    lx.kt[i,j]=sum(KT3[i,]!=KT3[j,])
}


row3.kt1=DCG.distance(KT3,DCG3.kt1,by.row=FALSE, k=7,method="euclidean",replicate=FALSE,r=3)

dist_row_large=row3.kt1+lx.kt
temp= c(0.04,0.07,0.1,0.25,0.5,10)

Ens3.kt2_3=Eigen.plot(temp, selected.id <- c( 1,2,3,4),dist_row_large)
DCG3.kt2_3=DCGtree.plot(num.clusters.selected=c(1,2,3,3,4,5),"mam 1st row tree",Ens3.kt2_3,temp)

save(Ens3.kt2_3, DCG3.kt2_3,file="KT_DCG2017_1.RData")

