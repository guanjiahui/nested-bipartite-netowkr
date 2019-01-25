load("DCG.RData")
nb=1000
overall=bootstrap.energy(1,1,DCG.10,DCG.11,BEEH,nb)
E.boot=numeric(nb)
BEEH.boot=list()
for (i in 1:nb){
  a=Bootbinary(BEEH)
  E.boot[i]=a$Energy
  BEEH.boot[[i]]=a$Matrix
}
save(overall,E.boot,BEEH.boot, file="overall_BEEH.RData")
#dat.comp=data.frame(dens=c(E.boot[1:nb],E.permute),lines=rep(c("bootstrap","permutation"),c(nb,1000)))
#densityplot(~dens,data=dat.comp,groups = lines,plot.points = FALSE, ref = TRUE,lty=c(2,1),
           # auto.key=list(lines=TRUE),main="B bootstrap comp")
