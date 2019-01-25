CLLO <- read.delim("~/Dropbox/Research/Pedro Jordan/Data Txt/CLLO.txt", header=FALSE)
cllo=CLLO[2:nrow(CLLO),2:ncol(CLLO)]
m2=nrow(cllo)
n2=ncol(cllo)
cllo=data.matrix(cllo)
cllo=matrix(as.numeric(cllo),nrow=m2)
ee=heatmap(cllo)
cllo=cllo[ee$rowInd,ee$colInd]
########
ly.cllo=matrix(0,n2,n2)
for (i in 1:n2){
  for(j in 1:n2)
    ly.cllo[i,j]=sum(cllo[,i]!=cllo[,j])
}
#########
#Apply DCG Tree
temp=c(0.8,2,5,12)
Ens.cllo=Eigen.plot(temp, selected.id=c( 1,2,3,4),ly.cllo)
DCG.cllo1=DCGtree.plot(num.clusters.selected=c(2,3,4,5),"cllo-initial column tree",Ens.cllo,temp)
#####################

row1.cllo=DCG.distance(cllo,DCG.cllo1,by.row=FALSE, k=3,method="euclidean",replicate=FALSE,r=1)
temp= c(0.1,0.15,0.4,1,1.5)
Ens.cllo2=Eigen.plot(temp, selected.id <- c( 1,2,3,4),row1.cllo)
DCG.cllo2=DCGtree.plot(num.clusters.selected=c(1,2,2,3,4),"cllo 1st row tree",Ens.cllo2,temp)
DCG2.cllo=DCGtree.plot(num.clusters.selected=c(2,3,4,8),"cllo 1st row tree",Ens.cllo2,temp)
############


########
GetBipEnergy(cllo)
GetBipEnergy(cllo[DCG2.cllo$order,DCG.cllo1$order])
##########
heatmap.2(cllo,Rowv=as.dendrogram(DCG.cllo2),margin=c(3,3),main="CLLO",
          Colv=as.dendrogram(DCG.cllo1),trace="none",col=scaleyellowred)
heatmap.2(cllo,Rowv=FALSE,Colv=FALSE,trace="none",col=scaleyellowred)
GetBipEnergy(cllo)
GetBipEnergy(cllo[DCG.cllo2$order,DCG.cllo1$order])
orr=cbind(DCG.cllo1$order,cutree(DCG.cllo1,k=10))
aaa=cutree(DCG.cllo2,k=10)
DCG.cllo2$order=which(cutree(DCG.cllo2,k=10))

############
library(geiger)




