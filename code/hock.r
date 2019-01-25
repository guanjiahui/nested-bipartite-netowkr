HOCK <- read.delim("~/Dropbox/Research/Pedro Jordan/Data Txt/HOCK.txt", header=FALSE)
hock=HOCK[2:nrow(HOCK),2:ncol(HOCK)]
m=nrow(hock)
n=ncol(hock)
hock=as.matrix(hock)
hock=matrix(as.numeric(hock),nrow=m)
heatmap(hock)
colnames(hock)=1:n
################
c=heatmap.2(hock,trace="none",col=scaleyellowred)
hock=hock[c$rowInd,c$colInd]
#########
ly.hock=matrix(0,n,n)
for (i in 1:n){
  for(j in 1:n)
    ly.hock[i,j]=sum(hock[,i]!=hock[,j])
}
heatmap(ly.hock,main="hamming distance")
levelplot(ly.hock)
#######
lx.hock=matrix(0,m,m)
for (i in 1:m){
  for(j in 1:m)
    lx.hock[i,j]=sum(hock[i,]!=hock[j,])
}

temp=c(0.4,0.7,1.5)
Ens.hockLX=Eigen.plot(temp,selected.id<-c(1,2,3),lx.hock)
h1=DCGtree.plot(num.clusters.selected=c(1,3,5),"hock",Ens.hockLX,temp)
######
col2.hock=DCG.distance(hock,h1,by.row=TRUE, k=8,method="euclidean",replicate=FALSE,r=5)
temp=c(0.15,1.2,5)
Ens_3=Eigen.plot(temp, selected.id <- c( 1,2,3),col2.hock)
h2=DCGtree.plot(num.clusters.selected=c(1,2,3),"hock-2nd column tree",Ens_3,temp)
######

#Apply DCG Tree
temp=c(0.5,1.2,5,12)
Ens.hock=Eigen.plot(temp, selected.id=c( 1,2,3,4),ly.hock)
DCG1.hock=DCGtree.plot(num.clusters.selected=c(2,2,3,5),"hock-initial column tree",Ens.hock,temp)
#####################
row1.hock=DCG.distance(hock,DCG1.hock,by.row=FALSE, k=5,method="euclidean",replicate=FALSE,r=5)
temp= c(0.2,0.3,0.4,1.5)
Ens.hock2=Eigen.plot(temp, selected.id <- c( 1,2,3,4),row1.hock)
DCG2.hock=DCGtree.plot(num.clusters.selected=c(1,2,3,5),"hock 1st row tree",Ens.hock2,temp)
############



col2.hock=DCG.distance(hock,DCG2.hock,by.row=TRUE, k=4,method="euclidean",replicate=FALSE,r=1)
temp=c(0.15,1.2,5)
Ens_hock3=Eigen.plot(temp, selected.id <- c( 1,2,3),col2.hock)
DCG3.hock=DCGtree.plot(num.clusters.selected=c(1,2,4),"hock-2nd column tree",Ens_hock3,temp)
######
row2.hock=DCG.distance(hock,DCG3.hock,by.row=FALSE, k=5,method="euclidean",replicate=FALSE,r=1)
temp= c(0.2,0.4,1.0)
Ens.hock4=Eigen.plot(temp, selected.id <- c( 1,2,3),row2.hock)
DCG4.hock=DCGtree.plot(num.clusters.selected=c(1,2,4),"2nd row tree",Ens.hock4,temp)

########


#####################
row1.hock=DCG.distance(hock,DCG1.hock,by.row=FALSE, k=5,method="euclidean",replicate=FALSE,r=1)
r1.hock=row1.hock+lx.hock
temp= c(0.2,0.3,0.4,1.5,10)
Ens.hock2=Eigen.plot(temp, selected.id <- c( 1,2,3,4,5),r1.hock)

D2.hock=DCGtree.plot(num.clusters.selected=c(1,3,9,15,19),"hock 1st row tree",Ens.hock2,temp)
##


temp= c(0.4,0.8,1.0,1.5,2,10)
Ens.hock2=Eigen.plot(temp, selected.id <- c( 1,2,3,4,5,6),r1.hock)

D2.hock2=DCGtree.plot(num.clusters.selected=c(1,3,4,9,13,18),"hock 1st row tree",Ens.hock2,temp)
############

HOCKrow=as.dendrogram(D2.hock2)
temp=HOCKrow[[2]][[1]]
HOCKrow[[2]][[1]]=HOCKrow[[1]]
HOCKrow[[1]]=temp




temp2=HOCKrow[[2]][[2]][[2]][[2]][[2]]
HOCKrow[[2]][[2]][[2]][[2]][[2]]=HOCKrow[[2]][[2]][[2]][[2]][[1]]
HOCKrow[[2]][[2]][[2]][[2]][[1]]=temp2

temp7=HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]
HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]=HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]

HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=temp7

tempt1=HOCKrow[[1]][[1]]
tempt2=HOCKrow[[1]][[2]][[1]]
HOCKrow[[1]][[1]]=HOCKrow[[1]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
HOCKrow[[1]][[2]][[1]]=tempt1
HOCKrow[[1]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=tempt2

temp37=HOCKrow[[1]][[2]][[1]]
temp15=HOCKrow[[1]][[2]][[2]][[1]]
HOCKrow[[1]][[2]][[1]]=HOCKrow[[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]
HOCKrow[[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]=temp15
HOCKrow[[1]][[2]][[2]][[1]]=temp37

###############
#tempThird=HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
#HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=
  #HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]

#HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=tempThird

#temptwo=HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
#HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=
 # HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
#HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=temptwo
###########
tempThird=HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=
  HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
HOCKrow[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=tempThird
#####################
temp1=HOCKcol2[[4]]
HOCKcol2[[4]]=HOCKcol2[[5]]
HOCKcol2[[5]]=temp1

temp3=HOCKcol2[[3]]
HOCKcol2[[3]]=HOCKcol2[[4]]
HOCKcol2[[4]]=temp3
###########################################################################
#check                                                                    #
heatmap.2(hock,Rowv=HOCKrow,Colv=rev(HOCKcol2),margin=c(2,18),            #
          trace = "none",main="HOCK",col=scaleyellowred)                  #
###########################################################################

heatmap.2(hock,Rowv=as.dendrogram(D2.hock2),Colv=rev(HOCKcol2),margin=c(2,18),           #
          trace = "none",main="HOCK",col=scaleyellowred)                                  #
######################

heatmap.2(hock[order.dendrogram(HOCKrow),rev(order.dendrogram(HOCKcol2))],Rowv=NULL,Colv=NULL,margin=c(2,18),           #
          trace = "none",main="HOCK",col=scaleyellowred)                                  #
######################
oo=rowSums(hock[order.dendrogram(HOCKrow),rev(order.dendrogram(HOCKcol2))])

GetBipEnergy(hock[D2.hock$order,DCG1.hock$order])
GetBipEnergy(hock[D2.hock2$order,DCG1.hock$order])
GetBipEnergy(hock[DCG2.hock$order,DCG1.hock$order])
GetBipEnergy(hock[DCG4.hock$order,DCG1.hock$order])
GetBipEnergy(hock[order.dendrogram(HOCKrow),rev(order.dendrogram(HOCKcol2))])




########################################
###############
library(gplots)



GetBipEnergy(hock[DCG2.hock$order,DCG1.hock$order])
GetBipEnergy(hock)
GetBipEnergy(hock[h1$order,h2$order])
h2order=c(1,2,3,4,h2$order[-(1:4)])
neworder=c(29, 28, 27, 26, 25, 24, 23, 22, 21 ,20 ,19 ,16 , 9 ,10,18 ,17 ,15 ,14, 13,
           12 , 8, 11 ,  7,  5 , 6,  4 , 3,  1 , 2)
hockDCG1=as.dendrogram(DCG1.hock)
hocknew=reorder(hockDCG1,order=neworder)
heatmap.2(hock,Rowv=as.dendrogram(DCG2.hock),Colv=hocknew,trace="none",col=scaleyellowred)
heatmap.2(hock,Rowv=FALSE,Colv=FALSE,trace="none",col=scaleyellowred)
######
heatmap.2(hock,Rowv=as.dendrogram(DCG2.hock),Colv=as.dendrogram(DCG1.hock),margin=c(3,4),
          trace = "none",main="HOCK",col=scaleyellowred)
########
colbranches <- function(n, col)
{
  a <- attributes(n) # Find the attributes of current node
  # Color edges with requested color
  attr(n, "edgePar") <- c(a$edgePar, list(col=col, lwd=2))
  n # Don't forget to return the node!
}
######
h2.hock=as.dendrogram(DCG2.hock)
h1.hock=as.dendrogram(DCG1.hock)

h22=as.phylo(DCG2.hock)
h11=as.phylo(DCG1.hock)
plot(h22,type="fan")
heatmap(hock,Rowv=h2.hock,Colv=h1.hock,trace="none",col=scaleyellowred)

Phy.hock=read.tree(text="(((((((5,40,39,2),((11,54),75),(70,10,58),(71,62),(78,46),(63,50,53,38,77,29,66,48,74,72,8)), 
              (9,7,4)),((((31,20,17,14),55),(65,28),((47,23),(27,22,3),
                                                     ((67,42),((80,56,49,34,30,21,18,15),32),
                                                                        64),26),51),70)),((60,37),((61,57,41,35),((12,25),33)))),
           ((76,24),((36,6),((59,1),16)))),68);")
plot(Phy.hock,direction="downwards",main="HOCK Animal Phylogeny")
plot(Phy.hock,type="unrooted")
#################
par(mfrow=c(1,2))
plot(Phy.hock,type="fan")
plot(h22,type="fan",tip.color=hsv(cutree(DCG2.hock,k=4)/4,1,1,0.7))
par(mfrow=c(1,1))
##############
orderh2=c(5,40,39,2,11,54,75,70,10,58,71,62,78,46,63,50,53,38,77,29,66,48,74,72,8, 
     9,7,4,31,20,17,14,55,65,28,47,23,27,22,3,67,42,80,56,49,34,30,21,18,15,32,64,26,51,70,60,37,61,57,41,35,12,25,33,
  76,24,36,6,59,1,16,68)
Phy.hock$tip.label=orderh2

Phy.hock2=read.tree(text="(((((22,14,15),(24,23)),10),((((21,20,16),18,9,12),((1,4),2),8),(28,29,26,17,13,3)),(((19,11),5,7,29),25)),6);")
plot(Phy.hock2,direction="downwards",, label.offset=0.4)
orderh1=as.numeric(Phy.hock2$tip.label)
#############
#fan plot 
par(mfrow=c(1,2))
plot(Phy.hock2,type="fan")
plot(h11,type="fan",tip.color=hsv(cutree(DCG1.hock,k=4)/4,1,1,0.7))
par(mfrow=c(1,1))
#######
######






HOCKcol=as.dendrogram(DCG1.hock)
ho1=cut(HOCKcol,h=1)
HOCKcol2=merge(ho1$lower[[1]],ho1$lower[[2]],ho1$lower[[3]][[1]],ho1$lower[[3]][[2]][[1]],ho1$lower[[3]][[2]][[2]])


###########################################################################################
heatmap.2(hock,Rowv=as.dendrogram(DCG2.hock),Colv=rev(HOCKcol2),margin=c(2,18),           #
          trace = "none",main="HOCK",col=scaleyellowred)                                  #
###########################################################################################


###########################
heatmap.2(hock,Rowv=as.dendrogram(DCG2.hock),Colv=rev(as.dendrogram(DCG1.hock)),margin=c(2,18),
          trace = "none",main="HOCK",col=scaleyellowred)
########
heatmap.2(hock[orderh2,orderh1],Rowv=FALSE,Colv=FALSE,margin=c(6,3),
          trace="none",col=scaleyellowred,main="HOCK with phylogeny")
#####
heatmap.2(hock[orderh2,DCG1.hock$order],Rowv=FALSE,Colv=FALSE,margin=c(6,3),
          trace="none",col=scaleyellowred,main="animal phylogeny")
heatmap.2(hock[DCG2.hock$order,orderh1],Rowv=FALSE,Colv=FALSE,margin=c(6,3),
          trace="none",col=scaleyellowred,main="plant phylogeny")
########
E1=GetBipEnergy(hock[DCG2.hock$order,DCG1.hock$order])
Ean=GetBipEnergy(hock[orderh2,DCG1.hock$order])
Epl=GetBipEnergy(hock[DCG2.hock$order,orderh1])
Eboth=GetBipEnergy(hock[orderh2,orderh1])
plot(c(E1,Ean,Epl,Eboth))
#####
n.random=5000
hock.permute=array(0,c(m,n,n.random))
for (i in 1:n.random){
  row.permute=sample(m,replace=FALSE)
  col.permute=sample(n,replace=FALSE)
  hock.permute[,,i]=hock[row.permute,]
  hock.permute[,,i]=hock.permute[,col.permute,i]
}
E.permute=sapply(1:n.random, function(i) GetBipEnergy(hock.permute[,,i]))##Energy of Permutation Dist
###
plot(density(E.permute),main="HOCK",xlim=c(-7650,-6250))
abline(v=E1,col=2)
abline(v=Ean,col=3)
abline(v=Epl,col=4)
abline(v=Eboth,col=5)
legend("topright",c("Coup Gemo","animal","plant","both"),col=c(2,3,4,5),lty=c(1,1,1,1),cex=0.6)
#############################################



hock.row=as.dendrogram(DCG2.hock)
hock.col=as.dendrogram(DCG1.hock)

row.order.hock=order.dendrogram(hock.row)
col.order.hock=order.dendrogram(hock.col)

hock.cp=hock[row.order.hock,col.order.hock]

hock.order=hock[order(rowSums(hock)),order(colSums(hock),decreasing=TRUE)]
hock.order2=hock.cp[order(rowSums(hock.cp)),order(colSums(hock.cp),decreasing=TRUE)]


heatmap.2(hock.order2,Rowv=FALSE,Colv=FALSE,  main="ordered row and col sum",margin=c(2,18),
          col=scaleyellowred,trace="none")


#########
#revised cp 
hock.cp=hock[order.dendrogram(HOCKrow),order.dendrogram(HOCKcol2)]
###############################################
################
r2=list()
r2[[1]]=1:5
r2[[2]]=6:7
r2[[3]]=8:25
r2[[4]]=26:28
r2[[5]]=29:31
r2[[6]]=32:44
r2[[7]]=45:77
r2[[8]]=78:81
#########
c2=list()
c2[[1]]=1:4
c2[[2]]=5:12
c2[[3]]=13:15
c2[[4]]=16:29
#########################################

CPforHOCK=N_CG(hock.cp,r2,c2,34,4,hock.cp)



########################

#########################


#x1=Miller_regular(hock.cp,5000)$Matrix
x1=list()
x1[[1]]=hock.cp
for (i in 2:6000){
  x1[[i]]=CheckerBoard(x1[[i-1]])
}


hE11=sapply(1:iter, function(i) GetBipEnergy(x1[[i]]))
#########
#delete the first depedenet matrices

hock.energy11=hE11[1001:5000]

hNCG11B=sapply(1001:6000,function(i) N_CG(x1[[i]],r2,c2,34,4,hock.cp))
library(bipartite)
hT11=sapply(1001:iter, function(i) nestedtemp(x1[[i]])$statistic)

hNODF11=sapply(1001:iter, function(i) nestednodf(x1[[i]],order=FALSE)[3]$statistic[3])
hNplus11=sapply(1001:iter,function(i) Nplus(x1[[i]]))

###########

#sample 2x2
iter=5000
sub1=Miller_regular(hock.cp[1:8,1:4],5000)$Matrix
sub2=Miller_regular(hock.cp[1:8,5:29],5000)$Matrix  

sub3=Miller_regular(hock.cp[9:16,1:4],5000)$Matrix  
sub4=Miller_regular(hock.cp[9:16,5:29],5000)$Matrix

sub5=Miller_regular(hock.cp[17:81,1:4],5000)$Matrix  
sub6=Miller_regular(hock.cp[17:81,5:29],5000)$Matrix

hE22=c()
hT22=c()
hNODF22=c()
hNplus22=c()

hNCG22B=c()

for (l in 1:iter){
  Binarytotal=rbind(cbind(sub1[,,l],sub2[,,l]),cbind(sub3[,,l],sub4[,,l]),cbind(sub5[,,l],sub6[,,l]))
  
  hE22[l]=GetBipEnergy(Binarytotal)
  hNODF22[l]=nestednodf(Binarytotal,order=FALSE)[3]$statistic[3]
  hT22[l]=nestedtemp(Binarytotal)$statistic
  hNplus22[l]=Nplus( Binarytotal)
  hNCG22B[l]=N_CG(Binarytotal,r2,c2,34,4,hock.cp)
}

#hNCG22=c()
#for (l in 1:5000){
 # Binarytotal=rbind(cbind(sub1[,,l],sub2[,,l]),cbind(sub3[,,l],sub4[,,l]))
  #hNCG22[l]=N_CG(Binarytotal,r2,c2,8,4,hock.cp)
#}
  

###########################

#sapmling 4x4

subm=matrix(list(),nrow=length(layer1$lower),ncol=4)

for (i in 1:iter){
  subm[i,j]=Miller_regular(hock.cp[rl1[[i]],cl1[[j]]],5000)$Matrix
}



iter=5000
sub11=Miller_regular(hock.cp[1:5,1:4],5000)$Matrix
sub12=Miller_regular(hock.cp[1:5,5:12],5000)$Matrix  
sub13=Miller_regular(hock.cp[1:5,13:15],5000)$Matrix  
sub14=Miller_regular(hock.cp[1:5,16:29],5000)$Matrix  

sub21=Miller_regular(hock.cp[6:28,1:4],5000)$Matrix
sub22=Miller_regular(hock.cp[6:28,5:12],5000)$Matrix  
sub23=Miller_regular(hock.cp[6:28,13:15],5000)$Matrix  
sub24=Miller_regular(hock.cp[6:28,16:29],5000)$Matrix  

sub31=Miller_regular(hock.cp[29:31,1:4],5000)$Matrix
sub32=Miller_regular(hock.cp[29:31,5:12],5000)$Matrix  
sub33=Miller_regular(hock.cp[29:31,13:15],5000)$Matrix  
sub34=Miller_regular(hock.cp[29:31,16:29],5000)$Matrix  

sub41=Miller_regular(hock.cp[32:81,1:4],5000)$Matrix
sub42=Miller_regular(hock.cp[32:81,5:12],5000)$Matrix  
sub43=Miller_regular(hock.cp[32:81,13:15],5000)$Matrix  
sub44=Miller_regular(hock.cp[32:81,16:29],5000)$Matrix  



hE44=c()
hT44=c()
hNODF44=c()
hNplus44=c()

hNCG44B=c()

for (l in 1:iter){
  Binarytotal=rbind(cbind(sub11[,,l],sub12[,,l],sub13[,,l],sub14[,,l]),
                    cbind(sub21[,,l],sub22[,,l],sub23[,,l],sub24[,,l]),
                    cbind(sub31[,,l],sub32[,,l],sub33[,,l],sub34[,,l]),
                    cbind(sub41[,,l],sub42[,,l],sub43[,,l],sub44[,,l]))
  hE44[l]=GetBipEnergy(Binarytotal)
  hNODF44[l]=nestednodf(Binarytotal,order=FALSE)[3]$statistic[3]
  hT44[l]=nestedtemp(Binarytotal)$statistic
  hNplus44[l]=Nplus( Binarytotal)
  hNCG44B[l]=N_CG(Binarytotal,r2,c2,8,4,hock.cp)
}

hNCG44=c()
for (l in 1:5000){
  Binarytotal=rbind(cbind(sub11[,,l],sub12[,,l],sub13[,,l],sub14[,,l]),
                    cbind(sub21[,,l],sub22[,,l],sub23[,,l],sub24[,,l]),
                    cbind(sub31[,,l],sub32[,,l],sub33[,,l],sub34[,,l]),
                    cbind(sub41[,,l],sub42[,,l],sub43[,,l],sub44[,,l]))
  hNCG44[l]=N_CG(Binarytotal,r2,c2,8,4,hock.cp)
}

#####################

plot(density(hock.energy11),col="dodgerblue",main="hock Energy",lwd=3,ylim=c(0,0.01),xlim=c(-7900,-7000))
lines(density(hE22),col="magenta",lwd=3)
lines(density(hE44),col="darkgreen",lwd=3)
abline(v=GetBipEnergy(hock.cp),lty=2,col="red",lwd=3)
abline(v=GetBipEnergy(hock.dcg.hc),lty=2,col="grey",lwd=3)
legend("topright",c("11","22","44", "DM"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))
#############################

plot(density(hNCG44),col="darkgreen",main="hock N_CG",lwd=3)
lines(density(hNCG22B),col="magenta",lwd=3)
lines(density(hNCG11B),col="dodgerblue",lwd=3)
abline(v=N_CG(hock.cp,r2,c2,8,4,hock.cp),lty=2,col="red",lwd=3)
legend("topright",c("11","22","44", "DM"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))
#############################

plot(density(hNODF11),col="dodgerblue",main="hock NODF",lwd=3,xlim=c(18,24),ylim=c(0,0.8))
lines(density(hNODF22),col="magenta",lwd=3)
lines(density(hNODF44),col="darkgreen",lwd=3)
abline(v=nestednodf(hock.cp,order=FALSE)[3]$statistic[3],lty=2,col="red",lwd=3)
legend("topright",c("11","22","44", "DM"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))
###

plot(density(hT11),col="dodgerblue",main="hock Temperature",lwd=3,ylim=c(0,2.2),xlim=c(3,7.5))
lines(density(hT22),col="magenta",lwd=3)
lines(density(hT44),col="darkgreen",lwd=3)
abline(v=nestedtemp(hock.cp)$statistic,lty=2,col="red",lwd=3)
legend("topright",c("11","22","44", "DM"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))
###

plot(density(hNplus11),col="dodgerblue",main="hock N Plus",lwd=3,xlim=c(90,250),ylim=c(0,0.028))
lines(density(hNplus22),col="magenta",lwd=3)
lines(density(hNplus44),col="darkgreen",lwd=3)
abline(v=Nplus(hock.cp),lty=2,col="red",lwd=3)
legend("topright",c("11","22","44", "DM"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))
################################


################################
###########################################################################
#check                                                                    #
heatmap.2(hock,Rowv=as.dendrogram(hc_row),Colv=rev(HOCKcol2),margin=c(2,18),            #
          trace = "none",main="HOCK",col=scaleyellowred)                  #
###########################################################################


fine=cut(HOCKrow,h=0.4)
r2=list()
for(i in 1:length(fine$lower)){
  members=labels(fine$lower[[i]])
  if (i==1)
    r2[[i]]=1:length(members)
  else
    r2[[i]]=(max(r2[[i-1]])+1):(max(r2[[i-1]])+length(members))
}
#####









#####################################
#comparison with the hierarchical clustering
########################################

hc_row=hclust(as.dist(lx.hock))
hc_col=hclust(as.dist(ly.hock))

#################################################################
#hc tree
heatmap.2(hock,Rowv=as.dendrogram(hc_row),Colv=as.dendrogram(hc_col),margin=c(2,18),      #
          trace = "none",main="HOCK HC Tree",col=scaleyellowred)   
#################################################################

hock.hc=hock[hc_row$order,hc_col$order]
hock.dcg.hc=hock[hc_row$order,order.dendrogram(HOCKcol2)]


#sample 2x2
iter=5000
sub1=Miller_regular(hock.hc[1:12,1:5],5000)$Matrix
sub2=Miller_regular(hock.hc[1:12,6:29],5000)$Matrix  

sub3=Miller_regular(hock.hc[13:81,1:5],5000)$Matrix  
sub4=Miller_regular(hock.hc[13:81,6:29],5000)$Matrix

hE22.hc=c()


for (l in 1:iter){
  Binarytotal=rbind(cbind(sub1[,,l],sub2[,,l]),cbind(sub3[,,l],sub4[,,l]))
  
  hE22.hc[l]=GetBipEnergy(Binarytotal)
}


###########################

#sapmling 4x4

iter=5000
sub11=Miller_regular(hock.hc[1:12,1:4],5000)$Matrix
sub12=Miller_regular(hock.hc[1:12,5:12],5000)$Matrix  
sub13=Miller_regular(hock.hc[1:12,13:15],5000)$Matrix  
sub14=Miller_regular(hock.hc[1:12,16:29],5000)$Matrix  

sub21=Miller_regular(hock.hc[13:28,1:4],5000)$Matrix
sub22=Miller_regular(hock.hc[13:28,5:12],5000)$Matrix  
sub23=Miller_regular(hock.hc[13:28,13:15],5000)$Matrix  
sub24=Miller_regular(hock.hc[13:28,16:29],5000)$Matrix  

sub31=Miller_regular(hock.hc[29:58,1:4],5000)$Matrix
sub32=Miller_regular(hock.hc[29:58,5:12],5000)$Matrix  
sub33=Miller_regular(hock.hc[29:58,13:15],5000)$Matrix  
sub34=Miller_regular(hock.hc[29:58,16:29],5000)$Matrix  

sub41=Miller_regular(hock.hc[59:81,1:4],5000)$Matrix
sub42=Miller_regular(hock.hc[59:81,5:12],5000)$Matrix  
sub43=Miller_regular(hock.hc[59:81,13:15],5000)$Matrix  
sub44=Miller_regular(hock.hc[59:81,16:29],5000)$Matrix  



hE44.hc=c()


for (l in 1:iter){
  Binarytotal=rbind(cbind(sub11[,,l],sub12[,,l],sub13[,,l],sub14[,,l]),
                    cbind(sub21[,,l],sub22[,,l],sub23[,,l],sub24[,,l]),
                    cbind(sub31[,,l],sub32[,,l],sub33[,,l],sub34[,,l]),
                    cbind(sub41[,,l],sub42[,,l],sub43[,,l],sub44[,,l]))
  hE44.hc[l]=GetBipEnergy(Binarytotal)

}


#######################

plot(density(hock.energy11),col="dodgerblue",main="hock",lwd=3,ylim=c(0,0.015),xlim=c(-7700,-7000))
lines(density(hE22.hc),col="magenta",lwd=3)
lines(density(hE44.hc),col="darkgreen",lwd=3)
abline(v=GetBipEnergy(hock.cp)+20,lty=2,col="red",lwd=3)
legend("topright",c("11","22","44", "HC"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))
#############################





save(hE44.hc,hE22.hc,hock.energy11,hock.hc,hc_row,hc_col,E34.hc,E22.hc,E11,BEEH.hc,
     file = "hc_hock_beeh.RData")




hock.dcg.hc=hock[hc_row$order,order.dendrogram(HOCKcol2)]



#################
#redo the DCG 
#################
n3=dim(hock.dcg.hc)[1]
ly2.hock.dcg.hc=matrix(0,n3,n3)
for (i in 1:n3){
  for(j in 1:n3)
    ly2.hock.dcg.hc[i,j]=sum(hock.dcg.hc[i,]!=hock.dcg.hc[j,])
}

temp=c(0.2,2,6)
Ens.hock.dcg.hc.ly=Eigen.plot(temp, selected.id=c( 1,2,3),ly2.hock.dcg.hc)

DCG.hock_row=DCGtree.plot(num.clusters.selected=c(1,2,7),
                      "hock.dcg.hc-initial row tree",Ens.hock.dcg.hc.ly,temp)


GetBipEnergy(hock.dcg.hc[DCG.hock_row$order,DCG.hock_col$order])

#################
#repeat the procedure. Now is get column distance from col tree 
col.hock.dcg.hc=DCG.distance(hock.dcg.hc,DCG.hock_row,by.row=TRUE, k=8,method="euclidean",replicate=FALSE,r=1)
temp=c(0.4,1.5,6)
Ens_hock.dcg.hc.col=Eigen.plot(temp, selected.id <- c( 1,2,3),col.hock.dcg.hc)

DCG.hock_col=DCGtree.plot(num.clusters.selected=c(1,2,4),"hock.dcg.hc-2nd column tree",Ens_hock.dcg.hc.col,temp)
######

#repeat the procedure. Now is get column distance from col tree 
row.hock.dcg.hc=DCG.distance(hock.dcg.hc,DCG.hock_col,by.row=FALSE, k=5,method="euclidean",replicate=FALSE,r=1)
row_input=ly2.hock.dcg.hc+row.hock.dcg.hc

temp=c(0.2,0.5,1,2,6)
Ens_hock.row=Eigen.plot(temp, selected.id <- c( 1,2,3,4,5),row_input)

DCG.hock_row2=DCGtree.plot(num.clusters.selected=c(1,2,3,6,7),"hock.dcg.hc-2nd column tree",Ens_hock.row,temp)
######


GetBipEnergy(hock.dcg.hc[DCG.hock_row2$order,DCG.hock_col$order])

GetBipEnergy(hock[hc_row$order,DCG.hock_col$order])
#hc tree
heatmap.2(hock.dcg.hc,Rowv=as.dendrogram(DCG.hock_row2),Colv=as.dendrogram(DCG.hock_col),margin=c(2,18),      #
          trace = "none",main="HOCK HC Tree",col=scaleyellowred)   

heatmap.2(hock.dcg.hc,Rowv=as.dendrogram(DCG.hock_row),Colv=as.dendrogram(DCG.hock_col),margin=c(2,18),      #
          trace = "none",main="HOCK HC Tree",col=scaleyellowred)   
#################################################################


hock.new=hock.dcg.hc[DCG.hock_row$order,DCG.hock_col$order]

#####################






plot(density(hock.energy11),col="dodgerblue",main="hock Energy",lwd=3,ylim=c(0,0.015),xlim=c(-7750,-7000))
lines(density(hE22.new),col="magenta",lwd=3)
lines(density(hE54.new),col="darkgreen",lwd=3)
lines(density(hE22.hc),col="magenta",lwd=3,lty=2)
lines(density(hE44.hc),col="darkgreen",lwd=3,lty=2)
abline(v=GetBipEnergy(hock.cp),lty=2,col="red",lwd=3)
#abline(v=GetBipEnergy(hock.dcg.hc),lty=2,col="grey",lwd=3)
abline(v=GetBipEnergy(hock.new),col="red",lwd=3)
legend("topright",c("11","22 DCG","54 DCG", "22 HC","44 HC","DM", "HC"),
       col=c("dodgerblue","magenta","darkgreen", "magenta","darkgreen","red","red"),
       lty=c(1,1,1,2,2,1,2),lwd=rep(3,7))
#########################################
























save(hE22.new,hE54.new,file="HOCK_En_Distrb.RData")






save(hock.new,DCG.hock_row,DCG.hock_col,DCG.hock_row2,file="HOCKnew.RData")



#sample 2x2
iter=5000
sub1=Miller_regular(hock.new[1:21,1:4],5000)$Matrix
sub2=Miller_regular(hock.new[1:21,5:29],5000)$Matrix  

sub3=Miller_regular(hock.new[22:81,1:4],5000)$Matrix  
sub4=Miller_regular(hock.new[22:81,5:29],5000)$Matrix

hE22.new=c()


for (l in 1:iter){
  Binarytotal=rbind(cbind(sub1[,,l],sub2[,,l]),cbind(sub3[,,l],sub4[,,l]))
  
  hE22.new[l]=GetBipEnergy(Binarytotal)
}

###########

iter=5000
sub11=Miller_regular(hock.new[1:21,1:4],5000)$Matrix
sub12=Miller_regular(hock.new[1:21,5:9],5000)$Matrix  
sub13=Miller_regular(hock.new[1:21,10:27],5000)$Matrix  
sub14=Miller_regular(hock.new[1:21,28:29],5000)$Matrix  

sub21=Miller_regular(hock.new[22:39,1:4],5000)$Matrix
sub22=Miller_regular(hock.new[22:39,5:9],5000)$Matrix  
sub23=Miller_regular(hock.new[22:39,10:27],5000)$Matrix  
sub24=Miller_regular(hock.new[22:39,28:29],5000)$Matrix  

sub31=Miller_regular(hock.new[40:57,1:4],5000)$Matrix
sub32=Miller_regular(hock.new[40:57,5:9],5000)$Matrix  
sub33=Miller_regular(hock.new[40:57,10:27],5000)$Matrix  
sub34=Miller_regular(hock.new[40:57,28:29],5000)$Matrix  

sub41=Miller_regular(hock.new[58:69,1:4],5000)$Matrix
sub42=Miller_regular(hock.new[58:69,5:9],5000)$Matrix  
sub43=Miller_regular(hock.new[58:69,10:27],5000)$Matrix  
sub44=Miller_regular(hock.new[58:69,28:29],5000)$Matrix  

sub51=Miller_regular(hock.new[70:81,1:4],5000)$Matrix
sub52=Miller_regular(hock.new[70:81,5:9],5000)$Matrix  
sub53=Miller_regular(hock.new[70:81,10:27],5000)$Matrix  
sub54=Miller_regular(hock.new[70:81,28:29],5000)$Matrix  

hE54.new=c()


for (l in 1:iter){
  Binarytotal=rbind(cbind(sub11[,,l],sub12[,,l],sub13[,,l],sub14[,,l]),
                    cbind(sub21[,,l],sub22[,,l],sub23[,,l],sub24[,,l]),
                    cbind(sub31[,,l],sub32[,,l],sub33[,,l],sub34[,,l]),
                    cbind(sub41[,,l],sub42[,,l],sub43[,,l],sub44[,,l]),
                    cbind(sub51[,,l],sub52[,,l],sub53[,,l],sub54[,,l]))
  hE54.new[l]=GetBipEnergy(Binarytotal)
  
}

hE54.new


#############################

N_CG=function(Mat,r,c,b1,b2,mam2)
{
  
  Block=matrix(list(),nrow=b1,ncol=b2)
  Intensity=matrix(0,nrow=b1,ncol=b2)
  
  for (i in 1:b1){
    for (j in 1:b2){
      Block[[i,j]]=Mat[r[[i]],c[[j]]]
      Intensity[i,j]=Tao_intensity(Block[[i,j]])
    }
  }#end of for i 
  
  # print(Intensity)
  
  # Block2=matrix(list(),nrow=b1,ncol=b2)
  #Intensity2=matrix(0,nrow=b1,ncol=b2)
  # 
  #for (i in 1:b1){
  #  for (j in 1:b2){
  #    Block2[[i,j]]=mam2[r[[i]],c[[j]]]
  #    Intensity2[i,j]=Tao_intensity(Block2[[i,j]])
  #  }
  # }#end of for i 
  
  #Sign=sign(Intensity-Intensity2)
  # Sign[which(Sign==0,arr.ind=TRUE)]=1
  
  # Intensity=Sign*Intensity
  # print(Intensity)
  
  C=expand.grid(1:b2,1:b2)
  C=C[-which(C$Var1<C$Var2),]
  DMT=0
  
  #compute the column ones (j)
  for (i in 1:b1){
    for(k in 1:nrow(C)){
      j=C[k,1]
      j1=C[k,2]
      w=Tao_intensity(mam2[r[[i]],])
      temp=(Intensity[i,j]-Intensity[i,j1])*(j-j1)*w
      Max=max(abs(Intensity[i,1:max(j,j1)]))
      
      
      if(temp==0 & Max==0)
        temp2=0
      else if (Max==0 & temp !=0)
        temp2=1e+6  # make it almost infinitive
      else 
        temp2=temp/Max
      
      DMT=temp+DMT
      #cat("j",j,"j1",j1,"i",i, "temp",temp,"Max",Max,"temp2",temp2,"\n")
    }
  } # end of for loop
  
  D=expand.grid(1:b1,1:b1)
  D=D[-which(D$Var1>D$Var2),]
  DMT2=0
  #next compute the row ones (i)
  
  for (j in 1:b2){
    for(k in 1:nrow(D)){
      i=D[k,1]
      i1=D[k,2]
      w=Tao_intensity(mam2[,c[[j]]])
      temp=(Intensity[i,j]-Intensity[i1,j])*(i1-i)*w
      Max=max(abs(Intensity[1:max(i,i1),j]))
      
      
      if(temp==0 & Max==0)
        temp2=0
      else if (Max==0 &temp !=0)
        temp2=1e+6  # make it almost infinitive
      else 
        temp2=temp/Max
      
      DMT2=temp+DMT2
      # cat("i",i,"i1",i1,"j",j, "temp",temp,"Max",Max,"temp2",temp2,"\n")
    }
  } # end of for loop
  ######
  
  ##########
  
  #######################################################################################
  
  DMT4=0
   C1=matrix(0,34,2)
  S_C=cbind(c(rep(-1,8),rep(1,26)), c(rep(-1,11),rep(1,23)))
  
  for (i in 1:b1){
    for (j in 3:b2){
      for (k in 2:(j-1)){
        
        delta_simulated=Intensity[i,(k+1)]+Intensity[i,(k-1)]-2*Intensity[i,k]
        # delta_DM=Intensity2[i,(k+1)]+Intensity2[i,(k-1)]-2*Intensity2[i,k]
        #s1=sign(delta_simulated)
        # s_DM[i,(k-1)]=sign(delta_DM)
        si=S_C[i,(k-1)]
        
        w=Tao_intensity(mam2[r[[i]],])
        #C1[i,(k-1)]=s1
        
        
        temp=-1*j*(nrow(Mat)-i+1)*si*(delta_simulated)*w
        DMT4=temp+DMT4
        
        #cat("i",i,"k",k-1,"Sign_DM",s1,"\n")
      }
    }
  }#end of for i 
  
  #print(C1)
  
  
  
  DMT3=0
  
  S_R=cbind(c(rep(-1,32)), c(rep(1,6),rep(-1,26)), c(rep(1,8),rep(-1,24)), c(rep(1,12),rep(-1,20)))
            
  
    R1=matrix(0,32,4)
  #cat("for DMT4 \n")
  
  for (i in 3:b1){
    for (j in 1:b2){
      for (k in 2:(i-1)){
        delta_simulated=Intensity[(k+1),j]+Intensity[(k-1),j]-2*Intensity[k,j]
        #delta_DM=Intensity2[(k+1),j]+Intensity2[(k-1),j]-2*Intensity2[k,j]
        #s1=sign(delta_simulated)
        # s_DM[(k-1),j]=sign(delta_DM)
        si=S_R[(k-1),j]
        #R1[(k-1),j]=s1
        
        w=Tao_intensity(mam2[,c[[j]]])
        temp=-1*j*(ncol(Mat)-j+1)*si*(delta_simulated)*w
        
        DMT3=temp+DMT3
        #cat("j",j,"k",k-1,"Sign_DM",s1, "temp",temp,"\n")
      }
    }
  }#end of for i 
  
  #print(R1)
  
  ###############
  DM_index=DMT+DMT2+DMT3+DMT4
  #cat("DMT",DMT,"DMT2",DMT2,"DMT3",DMT3,"DMT4",DMT4,"\n")
  
  return(DM_index)
  
  
}#DM function()


