#plot the energy density
load("DCG.RData")
library(lattice)
n.b=500
d1=bootstrap.energy(2,2,DCG.10,DCG.11,BEEH,n.b)
d2=bootstrap.energy(3,4,DCG.10,DCG.11,BEEH,n.b)
d3=bootstrap.energy(3,6,DCG.10,DCG.11,BEEH,n.b)
#d4=bootstrap.energy(6,3,DCG.10,DCG.11,BEEH,n.b)
############
save(d1,d2,d3,file="both_phylogeny.RData")
dat <- data.frame(dens = c(d1,d2,d3), lines = rep(c("coarse", "median","fine"), each = n.b))
densityplot(~dens,data=dat,groups = lines,plot.points = FALSE, ref = TRUE,lty=c(2,1,5),
            auto.key=list(lines=TRUE),main="BEEH bootstrap")
dat1 <- data.frame(dens = c(d2,d3), lines = rep(c("coarse", "fine"), each = n.b))
densityplot(~dens,data=dat1,groups = lines,plot.points = FALSE, ref = TRUE,lty=c(2,1),
            auto.key=list(lines=TRUE),main="BEEH bootstrap")