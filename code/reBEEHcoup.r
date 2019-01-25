#Set up the intial choice for ly (column) using Hamming Distance
m3=dim(BEEH)[1]
n3=dim(BEEH)[2]
ly.BEEH=matrix(0,n3,n3)
for (i in 1:n3){
  for(j in 1:n3)
    ly.BEEH[i,j]=sum(BEEH[,i]!=BEEH[,j])
}
heatmap(ly.BEEH,main="hamming distance")
levelplot(ly.BEEH)
#Apply DCG Tree
temp=c(0.5,0.72,1.2,5)
Ens.BEEH=Eigen.plot(temp, selected.id=c( 1,2,3,4),ly.BEEH)
DCG.1=DCGtree.plot(num.clusters.selected=c(1,2,3,4),"BEEH-initial column tree",Ens.BEEH,temp)
#########################
#Construct the row (DCG distance) based on the column tree 
row1.BEEH=DCG.distance(BEEH,DCG.1,by.row=FALSE, k=4,method="euclidean",replicate=FALSE,r=1)
temp= c(0.4,0.6,1.5)
Ens.BEEH2=Eigen.plot(temp, selected.id <- c( 1,2,3),row1.BEEH)
DCG.2=DCGtree.plot(num.clusters.selected=c(1,2,4),"BEEH 1st row tree",Ens.BEEH2,temp)
#########################

#repeat the procedure. Now is get column distance from row tree 
col2.BEEH=DCG.distance(BEEH,DCG.2,by.row=TRUE, k=3,method="euclidean",replicate=FALSE,r=1)
temp=c(0.1,0.4,1)
Ens_BEEH3=Eigen.plot(temp, selected.id <- c( 1,2,3),col2.BEEH)
DCG.3=DCGtree.plot(num.clusters.selected=c(1,3,5),"BEEH-2nd column tree",Ens_BEEH3,temp)
######
temp=c(0.4,0.45,0.5,0.55,1)
Ens_BEEH3=Eigen.plot(temp, selected.id <- c( 1,2,3,4,5),col2.BEEH)
DCG.3=DCGtree.plot(num.clusters.selected=c(1,2,4),"BEEH-2nd column tree",Ens_BEEH3,temp)
######


row2.BEEH=DCG.distance(BEEH,DCG.3,by.row=FALSE, k=4,method="euclidean",replicate=FALSE,r=1)
temp= c(0.4,0.6,1,1.5)
Ens.BEEH4=Eigen.plot(temp, selected.id <- c( 1,2,3,4),row2.BEEH)
DCG.4=DCGtree.plot(num.clusters.selected=c(1,2,5),"2nd row tree",Ens.BEEH4,temp)
####################
