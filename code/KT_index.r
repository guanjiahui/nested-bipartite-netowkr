load("KT1by1.RData")
library(bipartite)
iter=10000
T1by1.KT=sapply(1:iter,function(i) nestedtemp(Matt[[i]])$statistic)
NODF1by1.KT=sapply(1:iter,function(i) nestednodf(Matt[[i]],order=FALSE)[3]$statistic[3])
Matt_temp=lapply(1:iter, function(i) t(Matt[[i]])[rev(1:ncol(Matt[[i]])),])
Nplus1by1.KT=sapply(1:iter,function(i) Nplus(Matt_temp[[i]]))

save(T1by1.KT,NODF1by1.KT,Nplus1by1.KT, file="KTindex.RData")