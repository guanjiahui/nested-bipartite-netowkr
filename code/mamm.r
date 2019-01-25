mamm <- read.delim("~/Dropbox/Research/Mammalian/mamm.txt")
mam=mamm[,2:ncol(mamm)]
m1=nrow(mam)
n1=ncol(mam)
rownames(mam)=mamm[,1]
mam=data.matrix(mam)
mam=matrix(as.numeric(mam),nrow=m1)
d=heatmap(mam)
mam=mam[d$rowInd,d$colInd]
library(lattice)
####################
ly.mam=matrix(0,n1,n1)
for (i in 1:n1){
  for(j in 1:n1)
    ly.mam[i,j]=sum(mam[,i]!=mam[,j])
}
heatmap(ly.mam,main="hamming distance")
levelplot(ly.mam)
############
#Apply DCG Tree
temp=c(1.2,5,15)
Ens.mam=Eigen.plot(temp, selected.id=c( 1,2,3),ly.mam)
DCG.1=DCGtree.plot(num.clusters.selected=c(1,2,4),"mam-initial column tree",Ens.mam,temp)
########
row1.mam=DCG.distance(mam,DCG.1,by.row=FALSE, k=5,method="euclidean",replicate=TRUE,r=3)
temp= c(0.36,0.58,1.0,4)
Ens.mam2=Eigen.plot(temp, selected.id <- c( 1,2,3,4),row1.mam)
DCG.2=DCGtree.plot(num.clusters.selected=c(1,2,4,5),"mam 1st row tree",Ens.mam2,temp)
############

col2.mam=DCG.distance(mam,DCG.2,by.row=TRUE, k=5,method="euclidean",replicate=TRUE,r=5)
temp=c(0.28,0.53,1.5,5)
Ens.mam3=Eigen.plot(temp, selected.id=c( 1,2,3,4),col2.mam)
DCG.3=DCGtree.plot(num.clusters.selected=c(1,2,4,5),"mam-initial column tree",Ens.mam3,temp)
###########

scaleyellowred <- colorRampPalette(c("yellow", "black"), space = "rgb")(2)



heatmap.2(data.matrix(mamm[,2:ncol(mamm)]),Rowv=FALSE,Colv=FALSE, margin=c(3,7),trace="none",col=scaleyellowred)

heatmap.2(mam,Rowv=as.dendrogram(DCG.2),margin=c(3,7),main="EOLZ",
          Colv=as.dendrogram(DCG.3),trace="none",col=scaleyellowred)

heatmap.2(mam,Rowv=tree5,Colv=tree3,
          margin=c(3,3),trace="none",col=scaleyellowred)

heatmap.2(mamPerfect,Rowv=FALSE,Colv=FALSE,margin=c(3,3),trace="none",col=scaleyellowred)
heatmap.2(mam4,Rowv=FALSE,Colv=FALSE,margin=c(3,3),trace="none",col=scaleyellowred)
heatmap.2(t(LOW_matrix[[40]]),Rowv=FALSE ,Colv=FALSE,trace="none",col=scaleyellowred,main="lowest energy3")
############
############
tree5=as.dendrogram(DCG.2)
temp6=tree5[[2]][[2]][[2]]
tree5[[2]][[2]][[2]]=tree5[[2]][[2]][[1]]
tree5[[2]][[2]][[1]]=temp6

temp7=tree5[[2]][[1]][[2]]
tree5[[2]][[1]][[2]]=tree5[[2]][[1]][[1]]
tree5[[2]][[1]][[1]]=temp7

temp8=tree5[[2]][[1]][[1]][[2]][[1]]
tree5[[2]][[1]][[1]][[2]][[1]]=tree5[[2]][[1]][[1]][[2]][[2]][[2]]
tree5[[2]][[1]][[1]][[2]][[2]][[2]]=temp8

temp9=tree5[[2]][[1]][[1]][[2]][[1]]
tree5[[2]][[1]][[1]][[2]][[1]]=tree5[[2]][[1]][[1]][[2]][[2]][[1]]
tree5[[2]][[1]][[1]][[2]][[2]][[1]]=temp9
#check energy
library(gplots)
d=heatmap(mam)
GetBipEnergy(mam[d$rowInd,d$colInd])
GetBipEnergy(mam[DCG.2$order,DCG.3$order])
GetBipEnergy(mam[DCG.2$order,DCG.1$order])
GetBipEnergy(mam)
############
#permute the tree to make it looks as nested as possible
tree1=as.dendrogram(DCG.1)
flip_dend <- tree1;
tmp <- flip_dend[[2]][[1]]  # placeholder - 1,2 means "left, then right".
flip_dend[[2]][[1]] <- flip_dend[[2]][[2]][[2]][[1]] # move left branch to right
flip_dend[[2]][[2]][[2]][[1]] <- tmp; # move right branch to left.
##########
tree2=as.dendrogram(DCG.2)
temp=tree2[[2]][[2]][[2]][[1]]
tree2[[2]][[2]][[2]][[1]]=tree2[[2]][[2]][[2]][[2]]
tree2[[2]][[2]][[2]][[2]]=temp
#########
tree3=as.dendrogram(DCG.3)
temp=tree3[[2]][[2]][[2]][[2]]
tree3[[2]][[2]][[2]][[2]]=tree3[[2]][[1]]
tree3[[2]][[1]]=temp

temp2=tree3[[2]][[2]][[1]]
tree3[[2]][[2]][[1]]=tree3[[2]][[2]][[2]][[2]]
tree3[[2]][[2]][[2]][[2]]=temp2

temp3=tree3[[2]][[2]][[2]][[1]]
tree3[[2]][[2]][[2]][[1]]=tree3[[2]][[2]][[1]]
tree3[[2]][[2]][[1]]=temp3

temp4=tree3[[1]][[2]][[2]][[1]]
tree3[[1]][[2]][[2]][[1]]=tree3[[1]][[1]]
tree3[[1]][[1]]=temp4
###########
#bootstrap
row.order=c(3,5,4,1,2,9,7,8,25,24,23,6,22,28,27,26,21,20,19,18,17,16,14,15,13,12,10,11)
col.order=c(2,1,4,3,5,6,7,8,9,16,17,18,15,19,20,11,10,12,13,14,22,21,23,24,25,26)
col.order2=c(2,1,4,3,5,6,7,8,9,16,17,11,10,12,13,14,19,18,15,20,22,21,23,24,25,26)
library(lattice)
levelplot(t(t(mam)[row.order,col.order]))
mam2=t(mam)[row.order,col.order]
mam3=t(mam)[row.order,col.order2]
####
iter=200
Energy.coarse=numeric(iter)
Energy.fine=numeric(iter)
for (l in 1:iter){
  sub1=Bootbinary(mam2[1:5,1:20])$Matrix
  sub2=Bootbinary(mam2[1:5,21:26])$Matrix  
  
  sub3=Bootbinary(mam2[6:28,1:20])$Matrix  
  sub4=Bootbinary(mam2[6:28,21:26])$Matrix
  
  Binarytotal=rbind(cbind(sub1,sub2),cbind(sub3,sub4))
  Energy.coarse[l]=GetBipEnergy(Binarytotal)
}
######
dat.mam<- data.frame(dens = c(Emam22,Emam23,Emam551), lines = rep(c("E22","E23","E55"), c(200,200,500)))
densityplot(~dens,data=dat.mam,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammal")
####################
dat.mam2<- data.frame(dens = c(Emam22,Emam23,Emam25,Emam35,Emam52,Emam551,EProb11), 
                     lines = rep(c("E22","E23","E25","E35","E52","E55","E11"), 
                                 c(200,200,200,200,200,500,801)))

densityplot(~dens,data=dat.mam2,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammal")
#########
#perfec nestedness 
mamPN <- read.delim("~/Dropbox/Research/Mammalian/mamPN.txt", header=FALSE)
library(lattice)
mamPerfect=data.matrix(mamPN)
rownames(mamPerfect)=1:26
colnames(mamPerfect)=1:28
#############
#check whether we should merge two blocks 
A=mam3[,12:16]
B=mam3[,17:20]
AB=mam3[,12:20]
iter=200
EAB=numeric(iter)
for (i in 1:iter)
  EAB[i]=GetBipEnergy(ABmatrix[[i]])
#######
#1) Check the energy 
dat.AB<- data.frame(dens = c(E.AandB,EAB), lines = rep(c("A' and B'","A''B''"), c(200,200)))
densityplot(~dens,data=dat.AB,groups = lines,plot.points = FALSE, 
            ref =TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammal")
####################
#2) check the column sums
ColSumA=lapply(1:iter, function(i) colSums(t(ABmatrix[[i]][,1:5])))
ColSumAp=colSums(t(A))
ColSumB=lapply(1:iter, function(i) colSums(t(ABmatrix[[i]][,6:9])))
ColSumBp=colSums(t(B))
#####
ColSumA1=lapply(1:iter, function(i) ColSumA[[i]]+1)
ColSumAp1=ColSumAp+1
ColSumB1=lapply(1:iter, function(i) ColSumB[[i]]+1)
ColSumBp1=ColSumBp+1

ChisqrA=numeric(iter)
ChisqrB=numeric(iter)
for (i in 1:iter){
  ChisqrA[i]=sum((ColSumA1[[i]]-ColSumAp1)^2)
  ChisqrB[i]=sum((ColSumB1[[i]]-ColSumBp1)^2)
}
simulated=rchisq(200,df=25)
dat.chi<- data.frame(dens = c(ChisqrA,ChisqrB,simulated), lines = rep(c("Chisq A","Chisq B","Chisq simulated"), each=200))
densityplot(~dens,data=dat.chi,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="Chi-square")
####################
dat.AB2<- data.frame(dens = c(E.AandB2,EAB2), lines = rep(c("A' and B'","A''B''"), c(300,300)))
densityplot(~dens,data=dat.AB2,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammalSub")
####################
lower.limit <- min(dat.AB$dens)
upper.limit <- max(dat.AB$dens)
long.density <- density(subset(dat.AB, lines =="A' and B'")$dens, from = lower.limit, to = upper.limit, n = 2^10)
not.long.density <- density(subset(dat.AB, lines== "A''B''")$dens, from = lower.limit, to = upper.limit, n = 2^10)

density.difference <- long.density$y - not.long.density$y
intersection.point <- long.density$x[which(diff(density.difference > 0) != 0) + 1]

ggplot(dat.AB, aes(dens, fill = lines)) + geom_density(alpha = 0.2) + 
  geom_vline(xintercept = intersection.point, color = "red")

####

Total <- trapz (long.density$x, long.density$y) + trapz (long.density$x, not.long.density$y )
Overlap<-pmin(long.density$y,not.long.density$y)
Surface <- trapz ( long.density$x, Overlap ) / Total

##########
#for the second AB
lower.limit <- min(dat.AB2$dens)
upper.limit <- max(dat.AB2$dens)
long.density2 <- density(subset(dat.AB2, lines =="A' and B'")$dens, from = lower.limit, to = upper.limit, n = 2^10)
not.long.density2 <- density(subset(dat.AB2, lines== "A''B''")$dens, from = lower.limit, to = upper.limit, n = 2^10)

density.difference2<- long.density2$y - not.long.density2$y
intersection.point2 <- long.density2$x[which(diff(density.difference2 > 0) != 0) + 1]

ggplot(dat.AB2, aes(dens, fill = lines)) + geom_density(alpha = 0.2) + 
  geom_vline(xintercept = intersection.point2, color = "red")

####

Total <- trapz (long.density2$x, long.density2$y) + trapz (long.density2$x, not.long.density2$y )
Overlap<-pmin(long.density2$y,not.long.density2$y)
Surface2 <- trapz ( long.density2$x, Overlap) / Total
####################
#compute nestedness 
library(bipartite)
NEST=nestedness(mam3)
NEST2=nestedness(mam)
mamFull=matrix(1,26,28)
#####
iter=20
ETmam55=numeric(iter)
for (l in 109:200){
  ETmam55[l]=nestedness(Tmam55[[l]])$temperature 
}
#####
dat.T<- data.frame(dens = c(ETmam55,Tmam42,Tmam22,Tchange), lines = rep(c("T55","T42", "T22",""),each=200))
densityplot(~dens,data=dat.T,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammal T")
####################
T1mam55=list()
T1mam55=lapply(1:200, function(i) t(Tmam55[[i]])[rev(1:ncol(mam3)),])
NODF55=list()
NODF22=list()
for (l in 1:200){
  NODF55[[l]]=nestednodf(T1mam55[[l]],order=FALSE)
}
T2mam22=list()
T1mam22=lapply(1:200, function(i) t(Matrix22[[i]])[rev(1:ncol(mam3)),])
for (l in 1:200){
  NODF22[[l]]=nestednodf(T1mam22[[l]],order=FALSE)
}

T1mam42=lapply(1:200, function(i) t(Matrix42[[i]])[rev(1:ncol(mam3)),])
for (l in 1:200){
  NODF42[[l]]=nestednodf(T1mam42[[l]],order=FALSE)
}
nodf22=sapply(1:200, function(i) NODF22[[i]][3]$statistic[3])
nodf55=sapply(1:200, function(i) NODF55[[i]][3]$statistic[3])
nodf42=sapply(1:200, function(i) NODF42[[i]][3]$statistic[3])

dat.NODF<- data.frame(dens = c(nodf55,nodf42,nodf22), lines = rep(c("NODF55","NODF42", "NODF22"),each=200))
densityplot(~dens,data=dat.NODF,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammal NODF")
####################
Temp55Vegan=numeric(200)
for (l in 1:200)
  Temp55Vegan[l]=nestedtemp(Tmam55[[l]])$statistic
####
ind=which(ETmam55<3.7 & ETmam55>3.4)
ind.nodf=which(nodf55>47.5 & nodf55<48)
for (i in 1:length(ind))
heatmap.2(t(Tmam55[[47]]),Rowv=FALSE,Colv=FALSE,main="Both small bump #47",
          margin=c(3,3),trace="none",col=scaleyellowred)
heatmap.2(t(Tmam55[[93]]),Rowv=FALSE,Colv=FALSE,main="T small bump #93",
          margin=c(3,3),trace="none",col=scaleyellowred)
heatmap.2(t(Tmam55[[104]]),Rowv=FALSE,Colv=FALSE,main="NODF small bump #104",
          margin=c(3,3),trace="none",col=scaleyellowred)
############
heatmap.2(t(Tmam55[[50]]),Rowv=FALSE,Colv=FALSE,main="T majority #50",
          margin=c(3,3),trace="none",col=scaleyellowred)

 for (i in 1:length(index1))
A=sapply(1:length(index1), function(i) nestednodf(Tmam55[[index1[[i]]]],order=FALSE)[3]$statistic)

############
#permutation
n.random=200
mam.permute=array(0,c(nrow(mam2),ncol(mam2),n.random))
for (i in 1:n.random){
  row.permute=sample(nrow(mam2),replace=FALSE)
  col.permute=sample(ncol(mam2),replace=FALSE)
  mam.permute[,,i]=mam2[row.permute,]
  mam.permute[,,i]=mam.permute[,col.permute,i]
}
T.permute=sapply(1:n.random, function(i) nestednodf(mam.permute[,,i],order=FALSE)$statistic)
Tpermute=sapply(1:n.random,function(i) T.permute[,i][3])
heatmap.2(t(mam.permute[,,2]),Rowv=FALSE,Colv=FALSE,main="T majority",
          margin=c(3,3),trace="none",col=scaleyellowred)
dat.Tt<- data.frame(dens = c(Tmam55,Tmam42,Tmam22,Tpermute), lines = rep(c("T55","T42", "T22","Tpermute"),each=200))
densityplot(~dens,data=dat.Tt,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammal T")
####################
########
#2by2
T22=sapply(1:200, function(i) GetBipEnergy(Matrix22[[i]]))
heatmap.2((Matrix22[[38]]),Rowv=FALSE,Colv=FALSE,main="2by2 majority #38",
          margin=c(3,3),trace="none",col=scaleyellowred)

##############
#checker-board
Swap=function(Adj,i,j,i1,j1){
  swap=0
  if(Adj[i,j]==0 & Adj[i1,j]==1 & Adj[i,j1]==1 & Adj[i1,j1]==0){
    Adj[i,j]=1 
    Adj[i1,j]=0 
    Adj[i,j1]=0 
    Adj[i1,j1]=1
    swap=1
    return(list(Adj,swap))
  }
  if(Adj[i,j]==1 & Adj[i1,j]==0 & Adj[i,j1]==0 & Adj[i1,j1]==1){
    Adj[i,j]=0 
    Adj[i1,j]=1 
    Adj[i,j1]=1 
    Adj[i1,j1]=0
    swap=1
    return(list(Adj,swap))
  }# end if
  else
    return(list(Adj=Adj,swap=swap))
}#Swap
############
CheckerBoard=function(Adj,io,jo){
    for (i in 1: (nrow(Adj)-1)){
      for (j in 1:(ncol(Adj)-1)){
        if(i!=io & j!=jo){
        if (Swap(Adj,i,j)[[2]]==1){
          Matrix=Swap(Adj,i,j)[[1]]
          i1=i
          j1=j
          break
        }
        }
       }#end for j
    }#end for i
 # Energy=GetBipEnergy(Matrix)
  #Temperature=nestedtemp(Matrix)$statistic
 # NODF=nestednodf(Matrix,order=FALSE)[3]$statistic[3]
  return(list(Matrix=Matrix,i1=i1,j1=j1))#,Energy=Energy,Temperature=Temperature, NODF=NODF))
}
###########'

Checker=function(Adj,swap=1,iter,Mat=NULL){
  if(is.null(Mat))
    Mat<- list()
  
  if (swap==iter)
    return(Mat)
  #Stat=CheckerBoard(Adj,io,jo)
  mat=CheckerBoard(Adj)
  Mat[[swap]]=mat
  #i1=Stat$i1
  #j1=Stat$j1
  Mat<- Checker(mat,swap+1,iter,Mat)
  return(Mat)
}#Checker
####
iter=10000
BinaryChange=Checker(mam3,swap=1,iter)
Echange=sapply(1:(iter-1),function(i) GetBipEnergy(BinaryChange[[i]]))
NODFchange=sapply(1:(iter-1),function(i) nestednodf(BinaryChange[[i]],order=FALSE)[3]$statistic[3])
Tchange=sapply(1:(iter-1),function(i) nestedtemp(BinaryChange[[i]])$statistic)

###
CheckerBoard=function(Adj){
       i=1;j=1;i1=2;j1=2
        while(Swap(Adj,i,j,i1,j1)[[2]]!=1){
          A=sample(nrow(Adj),2,replace=TRUE)
          B=sample(ncol(Adj),2,replace=TRUE)
          i=A[1]
          j=B[1]
          i1=A[2]
          j1=B[2]
          Matrix=Swap(Adj,i,j,i1,j1)[[1]]
        }
  # Energy=GetBipEnergy(Matrix)
  #Temperature=nestedtemp(Matrix)$statistic
  # NODF=nestednodf(Matrix,order=FALSE)[3]$statistic[3]
  return(Matrix=Matrix)#,Energy=Energy,Temperature=Temperature, NODF=NODF))
}
####
iter=1500
BinaryChange2=Checker(mam3,swap=1,iter)
Echange2=sapply(1:(iter-1),function(i) GetBipEnergy(BinaryChange2[[i]]))
######
dat.NODF<- data.frame(dens = c(nodf55,nodf42,nodf22,NODFchange[200:1000]), lines = rep(c("NODF55","NODF42", "NODF22","NODFchekcer"),
                                                                                       c(200,200,200,801)))
densityplot(~dens,data=dat.NODF,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammal NODF")
####################
dat.T<- data.frame(dens = c(Tmam55,Tmam42,Tmam22,Tchange[200:1000]-1), lines = rep(c("T55","T42", "T22","Tchekcer"),
                                                                                 c(200,200,200,801)))
densityplot(~dens,data=dat.T,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="mammal T")
####################
plot(Echange,type="l")
abline(h=mean(Emam22),col="green")
abline(h=mean(Emam25),col="blue")
abline(h=mean(Emam55),col="red")
abline(h=mean(Echange[200:800]),col="orange")
legend("topright",c("55","25","22","11"),col=c("red","blue","green","orange"),lty=rep(1,4),cex=0.6)
########
#p = [26,24,23,21,19,13,13,12,11,10,10,9,9,9,7,7,7,7,6,6,5,5,4,3,2,1,1]
#q = [26,26,25,22,22,18,12,12,12,11,10,10,8,8,8,7,6,6,5,5,4,4,3,3,1,1]
########
#compare three different sampling (Patterson and Checker-Board)
dat.comp<- data.frame(dens = c(EPatterson[1:801],Echange[200:1000],Etemp[301:50000]), 
                      lines = rep(c("Patterson(800)","Checker(800)", "Checker(50000)"),
                               c(801,801, 49700)))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="Sampling Comparison")
####################
#N plus for 2x2, 4x2 and 5x5
Nplus22=sapply(1:200, function(i) Nplus((t(Matrix22[[i]])[rev(1:ncol(mam3)),])))
Nplus42=sapply(1:200, function(i) Nplus((t(Matrix42[[i]])[rev(1:ncol(mam3)),])))
Nplus55=sapply(1:200, function(i) Nplus((t(Tmam55[[i]])[rev(1:ncol(mam3)),])))


dat.comp<- data.frame(dens = c(Nplus22,Nplus42,Nplus55,Nchange[501:700]), 
                      lines = rep(c("22","42","55","11"), each=200))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Nplus",
            auto.key=list(lines=TRUE),main="Nplus")
####################

EE22=sapply(1:200, function(i) GetBipEnergy((t(Matrix22[[i]])[rev(1:ncol(mam3)),])))
EE42=sapply(1:200, function(i) GetBipEnergy((t(Matrix42[[i]])[rev(1:ncol(mam3)),])))
EE55=sapply(1:200, function(i) GetBipEnergy((t(Tmam55[[i]])[rev(1:ncol(mam3)),])))


dat.comp<- data.frame(dens = c(EE22,EE42,EE55), 
                      lines = rep(c("22","42","55"), each=200))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Nplus",
            auto.key=list(lines=TRUE),main="Nplus")
####################







##############################
#Use Miller's to do 1x1,2x2,4x2,5x5 for all the indexes (Energy, NODF, T, N+)
a <- rowSums(mam3)  # row sums
b <- colSums(mam3) # column sums
matrix_type <- 0  # 0: binary, 1: nonnegative integer

# Count the number of matrices with row sums a and column sums b
number_mam<- count1(a,b,matrix_type)

#####
number_of_samples <- 5000
x <- sample1(a,b,number_of_samples)
EMiller1x1=sapply(1:5000,function(i) GetBipEnergy(x[,,i]))
NODF.Miller1x1=sapply(1:5000,function(i) nestednodf(x[,,i],order=FALSE)[3]$statistic[3])
TMiller1x1=sapply(1:5000,function(i) nestedtemp(x[,,i])$statistic)
Nplus.Miller1x1=sapply(1:5000,function(i) Nplus( t(x[,,i])[rev(1:ncol(x[,,i])),]))
#################
#plot all the distribution 
#Energy
plot(EMiller1x1,type="l",ylab="Energy", main="Energy",
     col="blue",ylim=c(min(E.Miller.55)-80,max(EMiller1x1)))
#lines(EMiller1x1,col="grey")
#lines(Emam22,col="red")

#lines(E.Miller.42,col="green")
lines(E.Miller.55,col="red")
legend("bottomright",c("11 Miller","55"),col=c("blue","red"),lty=rep(1,2),cex=0.8)

########
#Energy Distribution 
dat.comp<- data.frame(dens = c(EMiller1x1,Emam22,E.Miller.42,E.Miller.55[201:400]), 
                      lines = rep(c("11","22","42","55"), c(rep(5000,3),200)))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="Mammal Energy")
####################

#Temperature Distribution 
dat.comp<- data.frame(dens = c(TMiller1x1,TMiller.22,TMiller.42,TMiller.55), 
                      lines = rep(c("11","22","42","55"), each=5000))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Temperature",
            auto.key=list(lines=TRUE),main="Mammal Temperature")
####################

#NODF Distribution 
dat.comp<- data.frame(dens = c(NODF.Miller1x1,NODF.Miller.22,NODF.Miller.42,NODF.Miller.55), 
                      lines = rep(c("11","22","42","55"), each=5000))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="NODF",
            auto.key=list(lines=TRUE),main="Mammal NODF")
####################

#N plus Distribution 
dat.comp<- data.frame(dens = c(Nplus.Miller1x1,Nplus.Miller.22,Nplus.Miller.42,Nplus.Miller.55[1:200]), 
                      lines = rep(c("11","22","42","55"), c(rep(5000,3),200)))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="N+",
            auto.key=list(lines=TRUE),main="Mammal N+")
################
par(mfrow=c(2,2))
plot(density(E.Miller.55),main="Energy 55")
plot(density(E.Miller.55[201:400]),main="Energy55, N=200")

plot(density(Nplus.Miller.55),main="N+ 11")
plot(density(Nplus.Miller.55[1:200]),main="N+ 11, N=200")
############
plot(density(E.Miller.55),main="Energy 55")
plot(density(E.Miller.55[201:400]),main="Energy55, N=200")

plot(density(EMiller1x1),main="Energy 11")
plot(density(EMiller1x1[1:200]),main="Energy11, N=200")

###########################################



########################
#Next index analysis
#make up the 5x5 blocks 
r=list()
r[[1]]=1:5
r[[2]]=6:8
r[[3]]=9:13
r[[4]]=14:24
r[[5]]=25:28
#########
c=list()
c[[1]]=1:2
c[[2]]=3:11
c[[3]]=12:15
c[[4]]=16:20
c[[5]]=21:26
#########

Block=matrix(list(),nrow=5,ncol=5)

for (i in 1:5){
  for (j in 1:5){
    Block[[i,j]]=mam2[r[[i]],c[[j]]]
  }
}
###############
DM1x1=sapply(1:5000,function(i) DM(t(x[,,i]),c,r,5,5))
DM21x1=sapply(1:5000,function(i) DM2(t(x[,,i]),c,r,5,5))
DM31x1=sapply(1:5000,function(i) DM3(t(x[,,i]),c,r,5,5))

DM41x1=sapply(1:5000,function(i) DM4(t(x[,,i]),c1,r1,5,5,t(mam3)))
EEEE=sapply(1:5000,function(i) GetBipEnergy(x[,,i]))
#### 
#DM distribution 
dat.comp<- data.frame(dens = c(DM1x1,DM.22,DM.42), 
                      lines = rep(c("11","22","42"), each=5000))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="G",
            auto.key=list(lines=TRUE),main="Mammal G")
#####################
dat.comp<- data.frame(dens = c(DM21x1,DM2.22,DM2.42), 
                      lines = rep(c("11","22","42"), each=5000))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="G",
            auto.key=list(lines=TRUE),main="Mammal G")
#########
pval1=(length(which(DM1x1<(-87.09747))))/5000
pval2=length(which(DM21x1<(-1625.536)))/5000

############################
dat.comp<- data.frame(dens = c(DM41x1,DM4.22,DM2.42), 
                      lines = rep(c("11","22","42"), each=5000))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="G",
            auto.key=list(lines=TRUE),main="Mammal G")

##################
GG11=sapply(1:5000,function(i) GG(t(x[,,i]),c1,r1,rowchoose,5,5))
dat.comp<- data.frame(dens = c(GG11,GG22,GG42), 
                      lines = rep(c("11","22","42"), each=5000))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="G",
            auto.key=list(lines=TRUE),main="Mammal G")
###########################################


ggplot(dat.comp, aes(dens, group = lines,colour=lines)) + geom_density(alpha=0.2) + 
  geom_vline(xintercept = 31, color = "red")

ggplot(dat.comp, aes(dens, group = lines,colour=lines)) + geom_density(alpha=0.2) + 
  geom_vline(xintercept = -691.2201, color = "red")

#####################################
#NCG index 

NCG11=sapply(1:5000,function(i) N_CG(t(x[,,i]),c1,r1,5,5,t(mam3)))

dat.comp<- data.frame(dens = c(NCG11,NCG22,NCG42), 
                      lines = rep(c("11","22","42"), each=5000))
densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="G",
            auto.key=list(lines=TRUE),main="Mammal G")
######
#refine plots 
# N_CG index 
plot(density(NCG42),xlim=c(min(NCG42)-20,max(NCG11)+40),col="darkgreen",xlab="N_CG",font.lab=2,font.axis=2,ylim=c(0,0.016),main="Mammal N_CG",lwd=3)
lines(density(NCG11),col="dodgerblue",lwd=3)
lines(density(NCG22),col="magenta",lwd=3)
abline(v=-399,lty=2,col="red",lwd=3)
legend("topright",c("11","22","42","DM"),col=c("dodgerblue","magenta","darkgreen","red"),lty=c(1,1,1,2),lwd=rep(3,4))
#############################

plot(density(E.Miller.55),xlim=c(min(E.Miller.55)-20,max(EMiller1x1)+20),col="orangered",xlab="Energy",ylim=c(0,0.028),font.lab=2,font.axis=2,main="Mammal Energy",lwd=3)
lines(density(EMiller1x1),col="dodgerblue",lwd=3)
lines(density(Emam22),col="magenta",lwd=3)
lines(density(E.Miller.42),col="darkgreen",lwd=3)
abline(v=-2184,lty=2,col="red",lwd=3)
legend("topright",c("11","22","42","55", "DM"),col=c("dodgerblue","magenta","darkgreen","orangered" ,"red"),lty=c(1,1,1,1,2),lwd=rep(3,5))
#############################

#NODF
plot(density(NODF.Miller.55),xlim=c(min(NODF.Miller1x1)-0.5,max(NODF.Miller1x1)+0.5),ylim=c(0,2.6),col="orangered",xlab="NODF",font.lab=2,font.axis=2,main="Mammal NODF",lwd=3)
lines(density(NODF.Miller1x1),col="dodgerblue",lwd=3)
lines(density(NODF.Miller.22),col="magenta",lwd=3)
lines(density(NODF.Miller.42),col="darkgreen",lwd=3)
abline(v=48.2062,lty=2,col="red",lwd=3)
legend("topleft",c("11","22","42","55", "DM"),col=c("dodgerblue","magenta","darkgreen","orangered" ,"red"),lty=c(1,1,1,1,2),lwd=rep(3,5))
#############################

#T
plot(density(TMiller.55-1),xlim=c(min(TMiller1x1-1)-0.5,max(TMiller1x1-1)+0.5),ylim=c(0,2.5),col="orangered",xlab="Temperature",font.lab=2,font.axis=2,main="Mammal Temperature",lwd=3)
lines(density(TMiller1x1-1),col="dodgerblue",lwd=3)
lines(density(TMiller.22-1),col="magenta",lwd=3)
lines(density(TMiller.42-1),col="darkgreen",lwd=3)
abline(v=3.581019,lty=2,col="red",lwd=3)
legend("topleft",c("11","22","42","55", "DM"),col=c("dodgerblue","magenta","darkgreen","orangered" ,"red"),lty=c(1,1,1,1,2),lwd=rep(3,5))
#############################

#N+
plot(density(Nplus.Miller.55),xlim=c(min(Nplus.Miller.22)-5,max(Nplus.Miller1x1)+5),ylim=c(0,0.22),col="orangered",xlab="N+",font.lab=2,font.axis=2,main="Mammal N+",lwd=3)
lines(density(Nplus.Miller1x1),col="dodgerblue",lwd=3)
lines(density(Nplus.Miller.22),col="magenta",lwd=3)
lines(density(Nplus.Miller.42),col="darkgreen",lwd=3)
abline(v=61,lty=2,col="red",lwd=3)
legend("topright",c("11","22","42","55", "DM"),col=c("dodgerblue","magenta","darkgreen","orangered" ,"red"),lty=c(1,1,1,1,2),lwd=rep(3,5))
#############################

plot(EMiller1x1,type="l",ylab="Energy", main="Energy",
     col="dodgerblue",ylim=c(min(E.Miller.55)-80,max(EMiller1x1)),font.lab=2,font.axis=2,xlab="Sampling index")
#lines(EMiller1x1,col="grey")
#lines(Emam22,col="red")

#lines(E.Miller.42,col="green")
lines(E.Miller.55,col="orangered")
legend("bottomright",c("11 Miller","55"),col=c("dodgerblue","orangered"),lty=rep(1,2),lwd=c(2,2))
##############################


#test for other co-occurrence indexes 
COMBO<-function(S){ncol(t(unique(t(S))))} 
############################
#1by1 indexes 
Checker_Diamond11=sapply(1:5000, function(i) checker(x[,,i]))
Cscore_Stone11=sapply(1:5000,function(i) c_score(x[,,i]))
Vratio_Roboson11=sapply(1:5000,function(i) V.ratio(x[,,i]))
Combo_Pielou11=sapply(1:5000,function(i) COMBO(x[,,i]))

#####
#other scale blocks indexes 
##########
#checker score index (Diamond 1975)

plot(density(Checker_Diamond22),xlim=c(min(Checker_Diamond11)-0.5,max(Checker_Diamond11)+0.5),col="magenta",xlab="checker(Diamond)",font.lab=2,font.axis=2,main="Mammal Checker(Diamond)",lwd=3)
lines(density(Checker_Diamond11),col="dodgerblue",lwd=3)
lines(density(Checker_Diamond42),col="darkgreen",lwd=3)
lines(density(Checker_Diamond55),col="orangered",lwd=3)
abline(v=6,lty=2,col="red",lwd=3)
legend("topright",c("11","22","42","55", "DM"),col=c("dodgerblue","magenta","darkgreen","orangered" ,"red"),lty=c(1,1,1,1,2),lwd=rep(3,5))
#############################

# C score (Stone and Roberts 1990)

plot(density(Cscore_Stone55),xlim=c(min(Cscore_Stone11)-0.2,max(Cscore_Stone11)+0.2),ylim=c(0,20),col="orangered",xlab="C score (Stone)",font.lab=2,font.axis=2,main="Mammal C score (Stone)",lwd=3)
lines(density(Cscore_Stone11),col="dodgerblue",lwd=3)
lines(density(Cscore_Stone22),col="magenta",lwd=3)
lines(density(Cscore_Stone42),col="darkgreen",lwd=3)
abline(v=2.661376,lty=2,col="red",lwd=3)
legend("topright",c("11","22","42","55", "DM"),col=c("dodgerblue","magenta","darkgreen","orangered" ,"red"),lty=c(1,1,1,1,2),lwd=rep(3,5))
#############################

# V ratio (Roboson 1972)

plot(Vratio_Roboson55,col="orangered",xlab="V ratio (Robson)",font.lab=2,font.axis=2,main="Mammal V ratio (Robson)",lwd=3)
lines(Vratio_Roboson11,col="dodgerblue",lwd=3)
lines(Vratio_Roboson22,col="magenta",lwd=3)
lines(Vratio_Roboson42,col="darkgreen",lwd=3)
abline(h=12.59493,lty=2,col="red",lwd=3)
legend("topright",c("11","22","42","55", "DM"),col=c("dodgerblue","magenta","darkgreen","orangered" ,"red"),lty=c(1,1,1,1,2),lwd=rep(3,5))
#############################

# COMBO  (Peilou 1969)
plot(density(Combo_Pielou55),xlim=c(min(Combo_Pielou42)-0.5,max(Combo_Pielou11)+0.5),ylim=c(0,7.5),col="orangered",xlab="COMBO (Pielou)",font.lab=2,font.axis=2,main="Mammal COMBO (Pielou)",lwd=3)
lines(density(Combo_Pielou11),col="dodgerblue",lwd=3)
lines(density(Combo_Pielou22),col="magenta",lwd=3)
lines(density(Combo_Pielou42),col="darkgreen",lwd=3)
abline(v=22,lty=2,col="red",lwd=3)
legend("topright",c("11","22","42","55", "DM"),col=c("dodgerblue","magenta","darkgreen","orangered" ,"red"),lty=c(1,1,1,1,2),lwd=rep(3,5))
#############################

sample(nrow(mam),nrow(mam))
m.random=mam[sample(nrow(mam),nrow(mam)),sample(ncol(mam),ncol(mam))]

#N+
plot(density(Nplus.Miller1x1),
     col="dodgerblue",xlab="Nestedness",font.lab=2,font.axis=2,
     main="Mammal Nestedness index",lwd=3)

abline(v=61,lty=2,col="red",lwd=3)
legend("topright",c("11","DM"),col=c("dodgerblue","red"),lty=c(1,2),lwd=rep(3,2))
#############################

plot(density(EMiller1x1),xlim=c(GetBipEnergy(mam3),max(EMiller1x1)),
     col="dodgerblue",xlab="Nestedness",font.lab=2,font.axis=2,
     main="Mammal Energy(Total Variation",lwd=3)

abline(v=GetBipEnergy(mam3),lty=2,col="red",lwd=3)
legend("topright",c("11","DM"),col=c("dodgerblue","red"),lty=c(1,2),lwd=rep(3,2))
#############################

heatmap.2(t(x[,,1]), trace="none",
          col=scaleyellowred,
          Rowv=FALSE,Colv=FALSE,main="sampled null matrix with fixed margin 1 ")

heatmap.2(t(x[,,10]), trace="none",
          col=scaleyellowred,
          Rowv=FALSE,Colv=FALSE,main="sampled null matrix with fixed margin 2")
