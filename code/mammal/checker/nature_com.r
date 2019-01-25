load("ecology.RData")
source("Nature_Communation.r")
source("Hfunction.r")

iter=1.0e+6
ran1=rnorm(length(mam2))
NatureEnsemble=lapply(1:iter,function(i) NatureComm(mam2))
HNature=sapply(1:iter, function(i) h1(NatureEnsemble[[i]],ran1))

#now check which two or three matrices are the same 
c1=unique(HNature)
ind.location=list()
tp=1
for (i in 1:length(c1)){
  ind.loc=which(HNature==c1[i])
  if (length(ind.loc)>1){
    ind.location[[tp]]=ind.loc
    tp=tp+1
  }
}

ind=unlist(ind.location)
MatRep_Nature=lapply(1:length(ind), function(i) NatureEnsemble[[ind[i]]])

save(HNature,ind.location,MatRep_Nature, file="Nature_comm.RData")
#######
