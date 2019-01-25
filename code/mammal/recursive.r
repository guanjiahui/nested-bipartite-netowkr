 load("checkerRevise.RData")
load("ecology.RData")
##########
c1=unique(Hv1)
ind.location=list()
tp=1
for (i in 1:length(c1)){
  ind.loc=which(Hv1==c1[i])
  if (length(ind.loc)>1){
    ind.location[[tp]]=ind.loc
    tp=tp+1
  }
  if (i==length(c1)/4)
    print("A quarter")
  if (i==length(c1)/2)
    print("half")
}

rep.time=sapply(1:length(ind.location),function(i) length(ind.location[[i]]))
xorder=sapply(1:length(ind.location),function(i) ind.location[[i]][1])
recursive_time=numeric(1e+6)
recursive_time[xorder]=1

save(ind.location, rep.time,xorder,recursive_time, file="recursive2.RData")