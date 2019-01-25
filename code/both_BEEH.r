#both phylogeny
load("DCG.RData")
library(lattice)
n.b=500
b1=bootstrap.energy(2,2,phy.an.tree,phy.pl.tree,BEEH,n.b)
b2=bootstrap.energy(3,4,phy.an.tree,phy.pl.tree,BEEH,n.b)
b3=bootstrap.energy(5,6,phy.an.tree,phy.pl.tree,BEEH,n.b)
save(b1,b2,b3,file="both_phylogeny1.Rdata")
dat.b=data.frame(dens=c(b1,b2),lines=rep(c("coarse","median","fine"),each=n.b))
densityplot(~dens,data=dat.b,groups = lines,plot.points = FALSE, ref = TRUE,lty=c(2,1,3),
            auto.key=list(lines=TRUE),main="BEEH bootstrap both col phylogeny")

#################################