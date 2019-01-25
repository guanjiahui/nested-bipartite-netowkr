EOLZ <- read.delim("~/Dropbox/Research/Pedro Jordan/Data Txt/EOLZ.txt", header=FALSE)
EZ=EOLZ[2:nrow(EOLZ),2:ncol(EOLZ)]
m1=nrow(EZ)
n1=ncol(EZ)
EZ=data.matrix(EZ)
EZ=matrix(as.numeric(EZ),nrow=m1)
d=heatmap(EZ)
EZ=EZ[d$rowInd,d$colInd]
####################
ly.EZ=matrix(0,n1,n1)
for (i in 1:n1){
  for(j in 1:n1)
    ly.EZ[i,j]=sum(EZ[,i]!=EZ[,j])
}
heatmap(ly.EZ,main="hamming distance")
levelplot(ly.EZ)
############
#Apply DCG Tree
temp=c(1.2,2.2,5,12)
Ens.EZ=Eigen.plot(temp, selected.id=c( 1,2,3,4),ly.EZ)
DCG.1=DCGtree.plot(num.clusters.selected=c(1,2,3),"EZ-initial column tree",Ens.EZ,temp)
DCG.11=DCGtree.plot(num.clusters.selected=c(2,2,3,4),"EZ-initial column tree",Ens.EZ,temp)
#####################
row1.EZ=DCG.distance(EZ,DCG.1,by.row=FALSE, k=3,method="euclidean",replicate=FALSE,r=1)
temp= c(0.12,0.22,0.8,2.5)
Ens.EZ2=Eigen.plot(temp, selected.id <- c( 1,2,3,4),row1.EZ)
DCG.2=DCGtree.plot(num.clusters.selected=c(1,2,3,5),"EZ 1st row tree",Ens.EZ2,temp)
############


Ens.EZ12=Eigen.plot(temp, selected.id <- c( 1,2,3,4),row1.EZ)
DCG.12=DCGtree.plot(num.clusters.selected=c(1,2,3,4),"EZ 1st row tree",Ens.EZ12,temp)

Ens.EZ22=Eigen.plot(temp, selected.id <- c( 1,2,3),row1.EZ)
DCG.22=DCGtree.plot(num.clusters.selected=c(2,2,3),"EZ 1st row tree",Ens.EZ22,temp)
Ens.EZ222=Eigen.plot(temp, selected.id <- c( 1,2,3,4),row1.EZ)
DCG.222=DCGtree.plot(num.clusters.selected=c(2,2,3,5),"EZ 1st row tree",Ens.EZ222,temp)

a=t(EOLZ)

b=EOLZ[,1]
colnames(EZ)=colnames(EOLZ)[-1]
rownames(EZ)=b[-1]
scaleyellowred <- colorRampPalette(c("black", "white"), space = "rgb")(2)
heatmap.2(t(EZ),Rowv=as.dendrogram(DCG.1),margin=c(3,3),main="EOLZ",
          Colv=as.dendrogram(DCG.12),trace="none",col=scaleyellowred)
heatmap.2(t(EZ),Rowv=FALSE,Colv=FALSE, margin=c(6,3),trace="none",col=scaleyellowred)
library(gplots)
d=heatmap(EZ)
GetBipEnergy(EZ)
GetBipEnergy(EZ[DCG.12$order,DCG.1$order])

GetBipEnergy(EZ[,DCG.1$order])
