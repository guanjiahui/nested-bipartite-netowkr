kt90 <- read.delim("~/Dropbox/Research/Pedro Jordan/phylogeny/kt90.txt", header=FALSE)
kt=kt90[2:nrow(kt90),2:ncol(kt90)]
m2=nrow(kt)
n2=ncol(kt)
kt=data.matrix(kt)
kt=matrix(as.numeric(kt),nrow=m2)
ee=heatmap(kt)
kt=kt[ee$rowInd,ee$colInd]
########
ly.kt=matrix(0,n2,n2)
for (i in 1:n2){
  for(j in 1:n2)
    ly.kt[i,j]=sum(kt[,i]!=kt[,j])
}
#########
#Apply DCG Tree
temp=c(0.8,2,7,15,80)
Ens.kt=Eigen.plot(temp, selected.id=c( 1,2,3,4,5),ly.kt)
DCG.kt1=DCGtree.plot(num.clusters.selected=c(1,2,3,4,6),"kt-initial column tree",Ens.kt,temp)
#####################

row1.kt=DCG.distance(kt,DCG.kt1,by.row=FALSE, k=4,method="euclidean",replicate=FALSE,r=1)
temp=c(0.8,2,7,15,80)
Ens.kt2=Eigen.plot(temp, selected.id <- c( 1,2,3,4),row1.kt)
DCG.kt2=DCGtree.plot(num.clusters.selected=c(1,1,1,2),"kt 1st row tree",Ens.kt2,temp)
############
heatmap.2(kt,Rowv=as.dendrogram(DCG3.kt2_2),Colv=as.dendrogram(DCG.kt1),margin=c(3,3),
          trace = "none",main="KT",col=scaleyellowred)
########
GetBipEnergy(kt)
GetBipEnergy(kt[DCG.kt2$order,DCG.kt1$order])
heatmap.2(kt,Rowv=FALSE,Colv=FALSE)

plot(Energy1by1.KT)
