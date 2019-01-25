CAFR.phy=CAFR_phylogeny[-1,-1]
CAFR.phy=as.matrix(CAFR.phy)
CAFR.phy=matrix(as.numeric(CAFR.phy),nrow=nrow(CAFR.phy))
m=nrow(CAFR.phy)
n=ncol(CAFR.phy)
for ( i in 1:m)
  for (j in 1:n)
    if (i<j)
    CAFR.phy[j,i]=CAFR.phy[i,j]
###########
CAFR.phy.plant=CAFR_phylogeny_plant[,-1]
CAFR.phy.plant=as.matrix(CAFR.phy.plant)
CAFR.phy.plant=matrix(as.numeric(CAFR.phy.plant),nrow=nrow(CAFR.phy.plant))
m2=nrow(CAFR.phy.plant)
n2=ncol(CAFR.phy.plant)
for ( i in 1:m2)
  for (j in 1:n2)
    if (i<j)
      CAFR.phy.plant[j,i]=CAFR.phy.plant[i,j]
#############
#apply DCG
animal=heatmap(CAFR.phy,main="CAFR.Phy animal")
plant=heatmap(CAFR.phy.plant,main="CAFR.phy plant")
phy.pl.CAFR=hclust(dist(CAFR.phy.plant))
phy.an.CAFR=hclust(dist(CAFR.phy))
phy.impact=Binary_CAFR[animal$rowInd,plant$rowInd]
levelplot(t(phy.impact),main="phylogeny impact")
heatmap(phy.impact,main="phylogeny impact")
########################
heatmap(phy.impact,Rowv=NA,Colv=NA)
heatmap(Binary_CAFR,Rowv=animal$rowInd,Colv=plant$rowInd)
#######

DCG11.dist=CAFR_DCG11[,-1]
DCG11.dist=as.matrix(DCG11.dist)
DCG11.dist=matrix(as.numeric(DCG11.dist),nrow=nrow(DCG11.dist))
m=nrow(DCG11.dist)
n=ncol(DCG11.dist)
for ( i in 1:m)
  for (j in 1:n)
    if (i<j)
      DCG11.dist[j,i]=DCG11.dist[i,j]
####

DCG12.dist=CAFR_DCG12[,-1]
DCG12.dist=as.matrix(DCG12.dist)
DCG12.dist=matrix(as.numeric(DCG12.dist),nrow=nrow(DCG12.dist))
m=nrow(DCG12.dist)
n=ncol(DCG12.dist)
for ( i in 1:m)
  for (j in 1:n)
    if (i<j)
      DCG12.dist[j,i]=DCG12.dist[i,j]