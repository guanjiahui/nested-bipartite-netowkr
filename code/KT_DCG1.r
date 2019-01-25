load("KT.RData")

row3.kt1=DCG.distance(KT3,DCG3.kt1,by.row=FALSE, k=5,method="euclidean",replicate=FALSE,r=1)
temp= c(0.2,0.8,10)
Ens3.kt2=Eigen.plot(temp, selected.id <- c( 1,2,3),row3.kt1)
DCG3.kt2=DCGtree.plot(num.clusters.selected=c(1,2,2),"mam 1st row tree",Ens3.kt2,temp)

save(DCG3.kt2,file="KT_DCG1.RData")