KT80 <- read.delim("~/Dropbox/Research/Mammalian/KT80/KT80.txt")
KT=KT80[,-1]
KT=data.matrix(KT)
A=rowSums(KT)
B=colSums(KT)
A1=order(A)
B1=order(B)
KT.order=KT[A1,B1]

heatmap.2(KT.order,trace="none")
heatmap.2(KT.order,Rowv=NULL,Colv=NULL,trace="none")

d=heatmap(KT)
KT=KT[d$rowInd,d$colInd]
######
m1=nrow(KT)
n1=ncol(KT)
ly.kt=matrix(0,n1,n1)
for (i in 1:n1){
  for(j in 1:n1)
    ly.kt[i,j]=sum(KT[,i]!=KT[,j])
}
levelplot(ly.kt)
############
#Apply DCG Tree
temp=c(0.5,1.2,2.5,6.5,15,60)
Ens.kt1=Eigen.plot(temp, selected.id=c( 1,2,3,4,5,6),ly.kt)
DCG.kt1=DCGtree.plot(num.clusters.selected=c(1,2,3,4,5,7),"mam-initial column tree",Ens.kt1,temp)
########
row1.kt=DCG.distance(KT,DCG.kt1,by.row=FALSE, k=4,method="euclidean",replicate=FALSE,r=1)
temp= c(0.5,1.0,10)
Ens.kt2=Eigen.plot(temp, selected.id <- c( 1,2,3),row1.kt)
DCG.kt2=DCGtree.plot(num.clusters.selected=c(1,2,2),"mam 1st row tree",Ens.kt2,temp)
############

#check 

scaleyellowred <- colorRampPalette(c("yellow", "black"), space = "rgb")(2)
heatmap.2(KT,Rowv=as.dendrogram(DCG.kt2),margin=c(3,7),main="KT",
          Colv=as.dendrogram(DCG.kt1),trace="none",col=scaleyellowred)
##########
d2=heatmap(KT.order)
KT2=KT.order[d2$rowInd,d2$colInd]
######
m1=nrow(KT)
n1=ncol(KT)
ly.kt2=matrix(0,n1,n1)
for (i in 1:n1){
  for(j in 1:n1)
    ly.kt2[i,j]=sum(KT2[,i]!=KT2[,j])
}

############
#Apply DCG Tree
temp=c(1.2,2.5,6.5,20,60)
Ens2.kt1=Eigen.plot(temp, selected.id=c( 1,2,3,4,5),ly.kt2)
DCG2.kt1=DCGtree.plot(num.clusters.selected=c(1,2,3,4,5),"mam-initial column tree",Ens2.kt1,temp)
########
row2.kt1=DCG.distance(KT2,DCG2.kt1,by.row=FALSE, k=4,method="euclidean",replicate=FALSE,r=1)
temp= c(0.2,0.8,10)
Ens2.kt2=Eigen.plot(temp, selected.id <- c( 1,2,3),row2.kt1)
DC2G.kt2=DCGtree.plot(num.clusters.selected=c(1,2,2),"mam 1st row tree",Ens2.kt2,temp)

DCG.kt22=DCGtree.plot(num.clusters.selected=c(1,2,4),"mam 1st row tree",Ens2.kt2,temp)
############
GetBipEnergy(KT2[DC2G.kt2$order,DCG2.kt1$order])
GetBipEnergy(KT2[DCG.kt22$order,DCG2.kt1$order])

heatmap.2(KT2,Rowv=as.dendrogram(DC2G.kt2),margin=c(3,7),main="KT",
          Colv=as.dendrogram(DCG2.kt1),trace="none",col=scaleyellowred)

heatmap.2(KT2,Rowv=as.dendrogram(DCG.kt22),margin=c(3,7),main="KT",
          Colv=as.dendrogram(DCG2.kt1),trace="none",col=scaleyellowred)

heatmap.2(KT2,Rowv=as.dendrogram(DC2G.kt2),margin=c(3,7),main="KT",
          Colv=TREE1,trace="none",col=scaleyellowred)

##########
TREE1=as.dendrogram(DCG2.kt1)
tmp1=TREE1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
TREE1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=TREE1[[1]][[1]][[1]]
TREE1[[1]][[1]][[1]]=tmp1

tmp2=TREE1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]
TREE1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]=TREE1[[1]][[1]][[2]][[1]]
TREE1[[1]][[1]][[2]][[1]]=tmp2


d2=heatmap(KT)
GetBipEnergy(KT[d2$rowInd,d2$colInd])
########
#start with hc of KT
KT3=KT[d2$rowInd,d2$colInd]
#############
m1=nrow(KT)
n1=ncol(KT)
ly.kt3=matrix(0,n1,n1)
for (i in 1:n1){
  for(j in 1:n1)
    ly.kt3[i,j]=sum(KT3[,i]!=KT3[,j])
}

#Apply DCG Tree
temp=c(1.2,2.5,6.5,20,80)
Ens3.kt1=Eigen.plot(temp, selected.id=c( 1,2,3,4,5),ly.kt3)
DCG3.kt1=DCGtree.plot(num.clusters.selected=c(2,2,3,4,5),"mam-initial column tree",Ens3.kt1,temp)
########
#check
load("~/Dropbox/Research/Bootstrap_server/KT/KT_DCG3.RData")
heatmap.2(KT3,Rowv=as.dendrogram(DCG3.kt2_3),margin=c(1,7),main="KT",
          Colv=as.dendrogram(DCG3.kt1),trace="none",col=scaleyellowred)
GetBipEnergy(KT3[DCG3.kt2_2$order,DCG3.kt1$order])
GetBipEnergy(KT3[DCG3.kt2_3$order,DCG3.kt1$order])
GetBipEnergy(KT3[DCG3.kt2_4$order,DCG3.kt1$order])
GetBipEnergy(KT3[DCG3.kt2_5$order,DCG3.kt1$order])
########
Tree.kt1=as.dendrogram(DCG3.kt2_3)
tmp1= Tree.kt1[[2]][[2]][[1]]
Tree.kt1[[2]][[2]][[1]]=Tree.kt1[[2]][[2]][[2]]
Tree.kt1[[2]][[2]][[2]]=tmp1

tmp2=Tree.kt1[[2]][[1]][[2]]
Tree.kt1[[2]][[1]][[2]]=Tree.kt1[[2]][[2]][[1]]
Tree.kt1[[2]][[2]][[1]]=tmp2

tmp3=Tree.kt1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]
Tree.kt1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]]=Tree.kt1[[1]][[1]][[1]]
Tree.kt1[[1]][[1]][[1]]=tmp3

tmp4=Tree.kt1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
Tree.kt1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]=Tree.kt1[[1]][[1]][[2]][[1]]
Tree.kt1[[1]][[1]][[2]][[1]]=tmp4

tmp5=Tree.kt1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]
Tree.kt1[[1]][[1]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[1]]= Tree.kt1[[1]][[1]][[2]][[2]][[1]]
Tree.kt1[[1]][[1]][[2]][[2]][[1]]=tmp5


GetBipEnergy(KT3[order.dendrogram(Tree.kt1),DCG3.kt1$order])

KT3.order1=KT3[order.dendrogram(Tree.kt1),DCG3.kt1$order]
KT3.order_original=KT3[DCG3.kt2_3$order,DCG3.kt1$order]

heatmap.2(KT,margin=c(1,7),main="KT",trace="none",col=scaleyellowred)
#################















############
#2016
temp= c(0.3,0.35,0.5,0.7,10,1000)
load("C:/Users/jiahui/Dropbox/Research/Bootstrap_server/KT/KT_DCG20163.RData")
D3.kt2=DCGtree.plot(num.clusters.selected=c(1,2,5,8,18,24),"mam 1st row tree",Ens3.kt2_2,temp)
########
KTRow=as.dendrogram(DCG3.kt2_2)

lay1=cut(KTRow,h=2.1)
d13=merge(lay1$lower[[1]],lay1$lower[[5]],lay1$lower[[6]],lay1$lower[[2]],lay1$lower[[4]],lay1$lower[[3]])
d13.r=Reduce(merge,d13)

d14=cut(d13,h=1.8)
d14part1=merge(d14$lower[[1]],d14$lower[[2]],d14$lower[[3]],d14$lower[[4]],
               d14$lower[[5]],d14$lower[[6]],d14$lower[[7]],d14$lower[[8]],
               d14$lower[[9]],d14$lower[[10]],d14$lower[[11]],d14$lower[[12]],
               d14$lower[[13]],d14$lower[[14]],d14$lower[[15]])
#d14Part1=Reduce(merge,d14part1)

d14p2=cut(d14$lower[[16]],h=1.2)
d14p22=merge(d14p2$lower[[2]],d14p2$lower[[4]],d14p2$lower[[3]],d14p2$lower[[1]])



d14part2=merge(d14$lower[[18]],d14$lower[[17]],d14p22)
d14Part2=Reduce(merge,d14part2)


###
dkt=merge(d14Part1,d14Part2)
KTrowTREE=Reduce(merge,dkt)
#check
ooo=rowSums(KT3[order.dendrogram(d13.r),DCG3.kt1$order])

heatmap.2(KT3,Rowv=as.dendrogram(DCG3.kt2_2),margin=c(2,7),main="KT",                                   #
          Colv=as.dendrogram(DCG3.kt1),trace="none",col=scaleyellowred)  

heatmap.2(KT3,Rowv=KTrowTREE,margin=c(2,7),main="KT",                                   #
          Colv=as.dendrogram(DCG3.kt1),trace="none",col=scaleyellowred)  


heatmap.2(KT3,Rowv=dkt,margin=c(2,7),main="KT",                                   #
          Colv=as.dendrogram(DCG3.kt1),trace="none",col=scaleyellowred)  

#######################################################################################
library(gplots)#
heatmap.2(KT3,Rowv=Tree.kt1,margin=c(2,7),main="KT",                                   #
          Colv=as.dendrogram(DCG3.kt1),trace="none",col=scaleyellowred)                #
                      

heatmap.2(KT3,Rowv=as.dendrogram(DCG3.kt2_3),margin=c(2,7),main="KT",                                   #
          Colv=as.dendrogram(DCG3.kt1),trace="none",col=scaleyellowred)                #
#
########################################################################################
KT.ordered=KT[order(rowSums(KT)),order(colSums((KT)),decreasing=TRUE)]
KT.ordered3=KT3[order(rowSums(KT3)),order(colSums((KT3)),decreasing=TRUE)]

heatmap.2(KT.ordered,Rowv=FALSE,margin=c(2,7),main="KT ordered row and column sums",
          Colv=FALSE,trace="none",col=scaleyellowred)
heatmap.2(KT.ordered3,Rowv=FALSE,margin=c(2,7),main="KT ordered row and column sums",
          Colv=FALSE,trace="none",col=scaleyellowred)

rtemp=rowSums(KT.ordered)
ctemp=colSums(KT.ordered)
r3temp=rowSums(KT.ordered3)
c3temp=colSums(KT.ordered3)
#################
levelplot(KT3.order1)
####

#1x1 mimicking using revised-checker-board
rand1=rnorm(length(KT3.order1))
iter=5000
#####
Hr2=c()
Hr2[1]=h1(KT3.order1,rand1)
ii=c();ii1=c()
jj=c();jj1=c()
Matt=list()
Matt[[1]]=KT3.order1
ii[1]=1;ii1[1]=1
jj[1]=2;jj1[1]=2
######
for (k in 2:iter){
  TMP=CheckerBoard_revision3(Matt[[k-1]],
                             ii[k-1],ii1[k-1],jj[k-1],jj1[k-1])
  
  Matt[[k]]=TMP$Matrix
  ii[k]=TMP$i
  ii1[k]=TMP$i1
  jj[k]=TMP$j
  jj1[k]=TMP$j1
  
  #Hr2[k]=h1(Matt[[k]],rand1)
}#end of for loop 
Energy1by1.KT=sapply(1:iter,function(i) GetBipEnergy(Matt[[i]]))
NODF1by1.KT=sapply(1:iter,function(i) nestednodf(Matt[[i]],order=FALSE)[3]$statistic[3])
T1by1.KT=sapply(1:iter,function(i) nestedtemp(Matt[[i]])$statistic)
Matt_temp=lapply(1:iter, function(i) t(Matt[[i]])[rev(1:ncol(Matt[[i]])),])
Nplus1by1.KT=sapply(1:iter,function(i) Nplus(Matt_temp[[i]]))

####################################################################################
#use Miller's by randomly selecting a rectangle matrix 

#first randomly choose the size of the selected matrix (length or rows and columns)
MillerRandom=function(Mat,iteration,n,a,b){
  EMiller=list()
  for (k in 1:iteration){
    u=0
    while(u==0){
      ro=sample(nrow(Mat),a)
      co=sample(ncol(Mat),b)
      
      
      matrix_type <- 0 
      KT.rand1=Mat[ro,co]
      a1=rowSums(KT.rand1)
      b1=colSums(KT.rand1)
      number_kt1 <- count1(a1,b1,matrix_type)
      number_kt1=as.numeric(number_kt1)
      
      if (number_kt1>1){
        u=1
        if (number_kt1<n){
          number_of_samples <- number_kt1
          KT1_sample<- sample1(a1,b1,number_of_samples)
        }
        else if (number_kt1>=n){
          number_of_samples <- n
          KT1_sample<- sample1(a1,b1,number_of_samples)
        }
      }#end of if that the count number of sample is not 1
      
    }#while loop
    
    KT.Miller=list()
    numSample=dim(KT1_sample)[3]
    for (i in 1:numSample){
      temp=Mat
      temp[ro,co]=KT1_sample[,,i]
      KT.Miller[[i]]=temp
    }
    
    
    EMiller[[k]]=sapply(1:number_of_samples,function(i) GetBipEnergy(KT.Miller[[i]]))
    
  }#end of for loop of total iteration of selecting "rectangles"
  return(EMiller)
}
###########################################################################
EnergyMiller=MillerRandom(Matt[[1000]],50,100,100,15)
EnergyMiller=unlist(EnergyMiller)
plot(c(Energy1by1.KT,EnergyMiller),type="l")

#############################################################################
#first randomly choose the size of the selected matrix (length or rows and columns)
MillerRandom2=function(Mat,n,nboc,INDEX,r=0,c=0,b1=0,b2=0,DM=NULL){

  KT_sample=list()
  
  ro=matrix(sample(nrow(Mat),nrow(Mat)),ncol=nboc)
  co=matrix(sample(ncol(Mat),ncol(Mat)),ncol=nboc)
  
  C=expand.grid(1:nboc,1:nboc)
  
  for (k in 1:nrow(C)){
    u=0
    while(u==0){

      matrix_type <- 0 
      KT.rand1=Mat[ro[,C[k,1]],co[,C[k,2]]]
      a1=rowSums(KT.rand1)
      b1=colSums(KT.rand1)
      number_kt1 <- count1(a1,b1,matrix_type)
      
      if (number_kt1>n){
        u=1

       number_of_samples <- number_kt1
       KT_sample[[k]]<- sample1(a1,b1,n)
  
      }#end of if that the count number of sample is not 1
      
    }#while loop
    
  }#end of for loop of total iteration of selecting "rectangles"
  
  KT.Miller=list()
  numSample=n
  for (i in 1:numSample){
    temp=Mat
    for (k in 1:nrow(C)){

      temp[ro[,C[k,1]],co[,C[k,2]]]=(KT_sample[[k]])[,,i]

    }
    KT.Miller[[i]]=temp
  }
  
  if (INDEX=="Energy")
    idx_value=sapply(1:n,function(i) GetBipEnergy(KT.Miller[[i]]))
  
  if (INDEX=="NCG")
    idx_value=sapply(1:n,function(i) N_CG(KT.Miller[[i]],r,c,b1,b2,DM))
  
  if (INDEX=="NODF")
    idx_value=sapply(1:n, function(i) nestednodf(KT.Miller[[i]],order=FALSE)[3]$statistic[3])
  
  if (INDEX=="Temperature")
    idx_value=sapply(1:n, function(i) nestedtemp(KT.Miller[[i]])$statistic)
  
  if (INDEX=="Nplus")
    idx_value=sapply(1:n,function(i) Nplus(KT.Miller[[i]]))
  else
    stop("please enter the correct index method")
  
  return(idx_value)
}


######################################################################







EMiller2=lapply(1:50,function(i) MillerRandom2(KT3.order1,100,6))
EMiller2=unlist(EMiller2)
plot(c(Energy1by1.KT,EMiller2),type="l")

EMiller2_3=lapply(1:50,function(i) MillerRandom2(KT3.order1,50,6))
EMiller2_3=unlist(EMiller2_3)


EMiller2_2=lapply(1:50,function(i) MillerRandom2(KT3.order1,50,5))
EMiller2_2=unlist(EMiller2_2)

##########################################

EMiller3=lapply(1:50,function(i) MillerRandom2(Matt[[1000]],100,6))
EMiller3=unlist(EMiller3)
EMiller3=EMiller3[-(3300:3400)]
plot(c(Energy1by1.KT,EMiller3),type="l")

dat1<- data.frame(dens = c(Energy1by1.KT[2001:10000],EMiller3,EMillerD[1001:5000]), 
                  lines = rep(c("Checker-Board","Independent Miller","Dependent Miller"), c(8000,length(EMiller3),4000)))
densityplot(~dens,data=dat1,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",lwd=3,font.lab=3,
            auto.key=list(lines=TRUE),main="Energy Distribution of sampled matrices")
####################
###############
plot(c(Energy1by1.KT,EMiller2_2,EMiller2_3),type="l",main="KT",ylab="Eneryg")
lines(10001:12500,EMiller2_2,col="blue")
lines(12501:15000,EMiller2_3,col="red")
legend("bottomright", c("Revised checker-board","Miller with 5 groups","Miller with 6 groups"),
       col=c("black","blue","red"),lty=c(1,1,1))

##########################












#########################
fine=cut(Tree.kt1,h=3)
r.kt=list()
for(i in 1:10){
  members=labels(fine$lower[[i]])
  if (i==1)
    r.kt[[i]]=1:length(members)
  else
    r.kt[[i]]=(max(r.kt[[i-1]])+1):(max(r.kt[[i-1]])+length(members))
}
#####
r.kt=list()
r.kt[[1]]=1:18
for (i in 2:9){
  members=labels(fine$lower[[i+1]])
  r.kt[[i]]=(max(r.kt[[i-1]])+1):(max(r.kt[[i-1]])+length(members))
}

###################
finec=cut(as.dendrogram(DCG3.kt1),h=0.1)
c.kt=list()
for (i in 1:7){
  members=labels(finec$lower[[i]])
  if (i==1)
    c.kt[[i]]=1:length(members)
  else
    c.kt[[i]]=(max(c.kt[[i-1]])+1):(max(c.kt[[i-1]])+length(members))
}

################

row.order.kt=order.dendrogram(Tree.kt1)
col.order.kt=order.dendrogram(as.dendrogram(DCG3.kt1))

kt.cp=KT3[row.order.kt,col.order.kt]




######################

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
  
  # cat("for DMG3 \n")
  DMT4=0
 # C1=matrix(0,9,5)
  S_C=cbind(c(-1,rep(1,8)),c(-1,-1,-1,rep(1,6)),c(rep(-1,5),rep(1,4)),
            c(rep(-1,7),1,1),c(rep(-1,8),1))
  for (i in 1:b1){
    for (j in 3:b2){
      for (k in 2:(j-1)){
        
        delta_simulated=Intensity[i,(k+1)]+Intensity[i,(k-1)]-2*Intensity[i,k]
        # delta_DM=Intensity2[i,(k+1)]+Intensity2[i,(k-1)]-2*Intensity2[i,k]
        s1=sign(delta_simulated)
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
  
  S_R=cbind(c(1,rep(-1,6)),c(1,1,-1,-1,-1,-1,-1),c(1,1,1,-1,-1,-1,-1),
            c(rep(1,4),rep(-1,3)),c(rep(1,5),-1,-1),c(rep(1,5),-1,-1), c(rep(1,6),-1))
#  R1=matrix(0,7,7)
  #cat("for DMT4 \n")
  for (i in 3:b1){
    for (j in 1:b2){
      for (k in 2:(i-1)){
        delta_simulated=Intensity[(k+1),j]+Intensity[(k-1),j]-2*Intensity[k,j]
        #delta_DM=Intensity2[(k+1),j]+Intensity2[(k-1),j]-2*Intensity2[k,j]
        s1=sign(delta_simulated)
        # s_DM[(k-1),j]=sign(delta_DM)
        #si=S_R[(k-1),j]
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

#########################
kt_NCGvalue=N_CG(kt.cp,r.kt,c.kt,9,7,kt.cp)

##############
#1by1 

KT11=MillerRandom2(Matt[[3000]],100,6,r.kt,c.kt,9,7,kt.cp)

KTNCG11=list()

for (i in 1:20){
   temp=MillerRandom2(Matt[[3000]],100,6,r.kt,c.kt,9,7,kt.cp)
   KTNCG11[[i]]=temp$NCG
   print(i)
}

##################################################

KNCG11=unlist(KTNCG11)

############################################################









#now doing 2x2
su1=list();su2=list();su3=list();su4=list()
su1[[1]]=kt.cp[1:42,1:12]
su2[[1]]=kt.cp[1:42,13:91]
su3[[1]]=kt.cp[43:679,1:12]
su4[[1]]=kt.cp[43:679,13:91]

for (i in 2:3000){
 # su1[[i]]=CheckerBoard(su1[[i-1]])
 # su2[[i]]=CheckerBoard(su2[[i-1]])
 # su3[[i]]=CheckerBoard(su3[[i-1]])
  su4[[i]]=CheckerBoard(su4[[i-1]])
}


kE22=c()
#kT22=c()
#kNODF22=c()
#kNplus22=c()

#kNCG22B=c()
iter=3000
for (l in 1500:iter){
  Binarytotal=rbind(cbind(su1[[l]],su2[[l]]),cbind(su3[[l]],su4[[l]]))
  kE22[l]=GetBipEnergy(Binarytotal)
 # kNODF22[l]=nestednodf(Binarytotal,order=FALSE)[3]$statistic[3]
 # kT22[l]=nestedtemp(Binarytotal)$statistic
 # kNplus22[l]=Nplus( Binarytotal)
 # kNCG22B[l]=N_CG(Binarytotal,r.kt,c.kt,9,7,kt.cp)
}

####################################################
#now doing 4x4
su11=list();su12=list();su13=list()
su21=list();su22=list();su23=list()
su31=list();su32=list();su33=list()
su41=list();su42=list();su43=list()

su14=list();su24=list();su34=list();su44=list()

su11[[1]]=kt.cp[1:17,1:12]
su12[[1]]=kt.cp[1:17,13:25]
su13[[1]]=kt.cp[1:17,26:46]
su14[[1]]=kt.cp[1:17,47:91]

su21[[1]]=kt.cp[18:42,1:12]
su22[[1]]=kt.cp[18:42,13:25]
su23[[1]]=kt.cp[18:42,26:46]
su24[[1]]=kt.cp[18:42,47:91]

su31[[1]]=kt.cp[43:87,1:12]
su32[[1]]=kt.cp[43:87,13:25]
su33[[1]]=kt.cp[43:87,26:46]
su34[[1]]=kt.cp[43:87,47:91]

su41[[1]]=kt.cp[88:679,1:12]
su42[[1]]=kt.cp[88:679,13:25]
su43[[1]]=kt.cp[88:679,26:46]
su44[[1]]=kt.cp[88:679,47:91]

for (i in 2:3000){
  
  su11[[i]]=CheckerBoard(su11[[i-1]])
  su12[[i]]=CheckerBoard(su12[[i-1]])
  su13[[i]]=CheckerBoard(su13[[i-1]])
  su14[[i]]=CheckerBoard(su14[[i-1]])
  
  su21[[i]]=CheckerBoard(su21[[i-1]])
  su22[[i]]=CheckerBoard(su22[[i-1]])
  su23[[i]]=CheckerBoard(su23[[i-1]])
  su24[[i]]=CheckerBoard(su24[[i-1]])
  
  su31[[i]]=CheckerBoard(su31[[i-1]])
  su32[[i]]=CheckerBoard(su32[[i-1]])
  su33[[i]]=CheckerBoard(su33[[i-1]])
  su34[[i]]=CheckerBoard(su34[[i-1]])
  
  su41[[i]]=CheckerBoard(su41[[i-1]])
  su42[[i]]=CheckerBoard(su42[[i-1]])
  su43[[i]]=CheckerBoard(su43[[i-1]])
  su44[[i]]=CheckerBoard(su44[[i-1]])
}

#######
kE44=c()
#kT44=c()
#kNODF44=c()
#kNplus44=c()

#kNCG44B=c()

for (l in 1500:iter){
  Binarytotal=rbind(cbind(sub11[[l]],sub12[[l]],sub13[[l]],sub14[[l]]),
                    cbind(sub21[[l]],sub22[[l]],sub23[[l]],sub24[[l]]),
                    cbind(sub31[[l]],sub32[[l]],sub33[[l]],sub34[[l]]),
                    cbind(sub41[[l]],sub42[[l]],sub43[[l]],sub44[[l]]))
  print(dim(Binarytotal))
  kE44[l]=GetBipEnergy(Binarytotal)
  
  
#  kNODF44[l]=nestednodf(Binarytotal,order=FALSE)[3]$statistic[3]
 # kT44[l]=nestedtemp(Binarytotal)$statistic
 # kNplus44[l]=Nplus( Binarytotal)
 # kNCG44B[l]=N_CG(Binarytotal,r.kt,c.kt,9,7,kt.cp)
}


for (l in 1500:iter){
  Binarytotal=rbind(cbind(su11[[l]],su12[[l]],su13[[l]],su14[[l]]),
                    cbind(su21[[l]],su22[[l]],su23[[l]],su24[[l]]),
                    cbind(su31[[l]],su32[[l]],su33[[l]],su34[[l]]),
                    cbind(su41[[l]],su42[[l]],su43[[l]],su44[[l]]))
  if(is.element(l, seq(1500,iter,by=100)))
  print(l)
  kE44[l]=GetBipEnergy(Binarytotal)
}




dat2<- data.frame(dens = c(EMiller3,kE22[1500:3000],kE44[1500:3000]), 
                  lines = rep(c("11","22","44"), c(length(EMiller3),length(kE22[1500:3000]),length(kE22[1500:3000]))))
densityplot(~dens,data=dat2,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",lwd=3,
            auto.key=list(lines=TRUE),main="Energy Distribution of KT")
####################













plot(c(Energy1by1.KT,EMiller3,EMillerD),type="l")

n_3=length(EMiller3)
plot(Energy1by1.KT[1:5000],col="dodgerblue",type="l",ylab="Energy",lwd=3,xlim=c(1,10000),font.lab=2,
     xlab = "Iteration",main="Eneregy of sampled matrices")
lines(EMillerD, col="magenta",lwd=3)
lines(1:n_3+5000,EMiller3, col="darkgreen",lwd=3)
legend("bottomright",c("checker-board","Dependent Miller","Indepedent Miller"), 
       lty=rep(1,3), col=c("dodgerblue","magenta","darkgreen"),lwd=3,cex=0.8)

#############



#Energy for KT
plot(density(EMiller3),col="dodgerblue",main="KT",lwd=3,ylim=c(0,0.005),xlab="Energy",
     xlim=c(-235650,-230200))
lines(density(kE22[1500:3000]),col="magenta",lwd=3)
lines(density(kE44[1500:3000]),col="darkgreen",lwd=3)
abline(v=GetBipEnergy(kt.cp),lty=2,col="red",lwd=3)
legend("topright",c("11","22","44", "DM"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))
#############################

#Energy for BEEH
plot(density(E34),col="darkgreen",main="BEEH",lwd=3,ylim=c(0,0.02),xlab="Energy",
     xlim=c(-600,-300))
lines(density(E11),col="dodgerblue",lwd=3)
lines(density(E22),col="magenta",lwd=3)
abline(v=-550,lty=2,col="red",lwd=3)
legend("topright",c("11","22","44", "DM"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))
#############################

#Energy for HOCK
plot(density(hock.energy11),col="dodgerblue",main="HOCK",lwd=3,xlab="Energy",
     ylim=c(0,0.01),xlim=c(-7700,-7000))
lines(density(hE22),col="magenta",lwd=3)
lines(density(hE44),col="darkgreen",lwd=3)
abline(v=GetBipEnergy(hock.cp),lty=2,col="red",lwd=3)
legend("topright",c("11","22","44", "DM"),col=c("dodgerblue","magenta","darkgreen", "red"),lty=c(1,1,1,2),lwd=rep(3,4))















#####################################
#comparison with the hierarchical clustering
########################################
m1=nrow(KT)
n1=ncol(KT)
lx.kt=matrix(0,m1,m1)
for (i in 1:m1){
  for(j in 1:m1)
    lx.kt[i,j]=sum(KT[i,]!=KT[j,])
}

hc_row_kt=hclust(as.dist(lx.kt))
hc_col_kt=hclust(as.dist(ly.kt))

#################################################################
#hc tree
heatmap.2(KT,Rowv=as.dendrogram(hc_row_kt),Colv=as.dendrogram(hc_col_kt),     #
          trace = "none",main="KT HC Tree",col=scaleyellowred)   
#################################################################
KT.hc=KT[hc_row_kt$order,hc_col_kt$order]


#now doing 2x2
su1=list();su2=list();su3=list();su4=list()
su1[[1]]=KT.hc[1:330,1:2]
su2[[1]]=KT.hc[1:330,3:91]
su3[[1]]=KT.hc[331:679,1:2]
su4[[1]]=KT.hc[331:679,3:91]

for (i in 1250:3000){
   #su1[[i]]=CheckerBoard(su1[[i-1]])
   #su2[[i]]=CheckerBoard(su2[[i-1]])
  # su3[[i]]=CheckerBoard(su3[[i-1]])
  su4[[i]]=CheckerBoard(su4[[i-1]])
}


kE22.hc=c()

iter=3000
for (l in 1500:iter){
  Binarytotal=rbind(cbind(su1[[l]],su2[[l]]),cbind(su3[[l]],su4[[l]]))
  kE22.hc[l]=GetBipEnergy(Binarytotal)
  # kNODF22[l]=nestednodf(Binarytotal,order=FALSE)[3]$statistic[3]
  # kT22[l]=nestedtemp(Binarytotal)$statistic
  # kNplus22[l]=Nplus( Binarytotal)
  # kNCG22B[l]=N_CG(Binarytotal,r.kt,c.kt,9,7,kt.cp)
}

save(kE22.hc,KT.hc,file="KTHC_server.RData")
######################################
#now doing 4x4
su11=list();su12=list();su13=list()
su21=list();su22=list();su23=list()
su31=list();su32=list();su33=list()
su41=list();su42=list();su43=list()

su14=list();su24=list();su34=list();su44=list()

su11[[1]]=KT.hc[1:30,1:2]
su12[[1]]=KT.hc[1:30,3:30]
su13[[1]]=KT.hc[1:30,31:80]
su14[[1]]=KT.hc[1:30,81:91]

su21[[1]]=KT.hc[31:330,1:2]
su22[[1]]=KT.hc[31:330,3:30]
su23[[1]]=KT.hc[31:330,31:80]
su24[[1]]=KT.hc[31:330,81:91]

su31[[1]]=KT.hc[331:650,1:2]
su32[[1]]=KT.hc[331:650,3:30]
su33[[1]]=KT.hc[331:650,31:80]
su34[[1]]=KT.hc[331:650,81:91]

su41[[1]]=KT.hc[651:679,1:2]
su42[[1]]=KT.hc[651:679,3:30]
su43[[1]]=KT.hc[651:679,31:80]
su44[[1]]=KT.hc[651:679,81:91]

for (i in 2:3000){
  
  su11[[i]]=CheckerBoard(su11[[i-1]])
  su12[[i]]=CheckerBoard(su12[[i-1]])
  su13[[i]]=CheckerBoard(su13[[i-1]])
  su14[[i]]=CheckerBoard(su14[[i-1]])
  
  su21[[i]]=CheckerBoard(su21[[i-1]])
  su22[[i]]=CheckerBoard(su22[[i-1]])
  su23[[i]]=CheckerBoard(su23[[i-1]])
  su24[[i]]=CheckerBoard(su24[[i-1]])
  
  su31[[i]]=CheckerBoard(su31[[i-1]])
  su32[[i]]=CheckerBoard(su32[[i-1]])
  su33[[i]]=CheckerBoard(su33[[i-1]])
  su34[[i]]=CheckerBoard(su34[[i-1]])
  
  su41[[i]]=CheckerBoard(su41[[i-1]])
  su42[[i]]=CheckerBoard(su42[[i-1]])
  su43[[i]]=CheckerBoard(su43[[i-1]])
  su44[[i]]=CheckerBoard(su44[[i-1]])
}

#######
kE44=c()



for (l in 1500:iter){
  Binarytotal=rbind(cbind(su11[[l]],su12[[l]],su13[[l]],su14[[l]]),
                    cbind(su21[[l]],su22[[l]],su23[[l]],su24[[l]]),
                    cbind(su31[[l]],su32[[l]],su33[[l]],su34[[l]]),
                    cbind(su41[[l]],su42[[l]],su43[[l]],su44[[l]]))
  if(is.element(l, seq(1500,iter,by=100)))
    print(l)
  kE44[l]=GetBipEnergy(Binarytotal)
}









################
CheckerBoard=function(Adj){
  i=1;j=1;i1=2;j1=2
  u=0
  while(u!=1){
    
    A=sample(nrow(Adj),2,replace=TRUE)
    B=sample(ncol(Adj),2,replace=TRUE)
    i=A[1]
    j=B[1]
    i1=A[2]
    j1=B[2]
    
    u=Swap(Adj,i,j,i1,j1)[[2]]
    matrix.result=Swap(Adj,i,j,i1,j1)[[1]]
  }

  return(matrix.result)
}
