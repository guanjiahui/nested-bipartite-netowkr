load("DCG.RData")
library(lattice)
n.b=500
c1=bootstrap.energy(2,2,phy.an.tree,DCG.1,BEEH,n.b)
c2=bootstrap.energy(3,4,phy.an.tree,DCG.1,BEEH,n.b)
c3=bootstrap.energy(3,7,phy.an.tree,DCG.1,BEEH,n.b)
########
save(c1,c2,c3,file="row_phylogeny.RData")
dat.c=data.frame(dens=c(c1,c2,c3),lines=rep(c("coarse","median","fine"),each=n.b))
densityplot(~dens,data=dat.c,groups = lines,plot.points = FALSE, ref = TRUE,lty=c(2,1,3),
            auto.key=list(lines=TRUE),main="BEEH bootstrap fixed row phylogeny")