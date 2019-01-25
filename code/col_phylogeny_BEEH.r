load("DCG.RData")
library(lattice)
n.b=500
r1=bootstrap.energy(2,2,DCG.1.lx,phy.pl.tree,BEEH,n.b)
r2=bootstrap.energy(3,4,DCG.1.lx,phy.pl.tree,BEEH,n.b)
save(r1,r2,file="row_phylogeny.RData")
#######
dat.r=data.frame(dens=c(r1,r2),lines=rep(c("coarse","fine"),each=n.b))
densityplot(~dens,data=dat.r,groups = lines,plot.points = FALSE, ref = TRUE,lty=c(2,1),
            auto.key=list(lines=TRUE),main="BEEH bootstrap fixed col phylogeny")