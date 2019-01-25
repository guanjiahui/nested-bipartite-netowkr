FROS.phy=FROS_PL.phylogeny[,-1]
FROS.phy=as.matrix(FROS.phy)
FROS.phy=matrix(as.numeric(FROS.phy),nrow=nrow(FROS.phy))
m=nrow(FROS.phy)
n=ncol(FROS.phy)
for ( i in 1:m)
  for (j in 1:n)
    if (i<j)
      FROS.phy[j,i]=FROS.phy[i,j]
###########
FROS.phy.plant=FROS_PL.phylogeny[-1,-1]
b=colnames(FROS.phy.plant)
FROS.phy.plant=as.matrix(FROS.phy.plant)
FROS.phy.plant=matrix(as.numeric(FROS.phy.plant),nrow=nrow(FROS.phy.plant))
m2=nrow(FROS.phy.plant)
n2=ncol(FROS.phy.plant)
for ( i in 1:m2)
  for (j in 1:n2)
    if (i<j)
      FROS.phy.plant[j,i]=FROS.phy.plant[i,j]
rownames(FROS.phy.plant)=b
rownames(FROS.phy)=rownames(FROS)
############
################
animal=heatmap(FROS.phy,main="FROS.Phy animal",margin=c(5,11))
plant=heatmap(FROS.phy.plant,main="FROS.phy animal")
phy.pl.FROS=hclust(dist(FROS.phy.plant))