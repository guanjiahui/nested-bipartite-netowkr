load("ecology.RData")
source("Hfunction.r")
source("checkerBoard.r")

###########
rand1=rnorm(length(mam2))
iter=1.0e+6
#####
Hr3=c()
Hr3[1]=h1(mam2,rand1)
ii=c();ii1=c()
jj=c();jj1=c()
Matt=list()
Matt[[1]]=mam2
ii[1]=1;ii1[1]=1
jj[1]=2;jj1[1]=2
######
for (k in 2:iter){
  TMP=CheckerBoard_revision4(Matt[[k-1]],
                             ii[k-1],ii1[k-1],jj[k-1],jj1[k-1])
  
  Matt[[k]]=TMP$Matrix
  ii[k]=TMP$i
  ii1[k]=TMP$i1
  jj[k]=TMP$j
  jj1[k]=TMP$j1
  
  Hr3[k]=h1(Matt[[k]],rand1)
}#end of for loop 

###########
#now check which two or three matrices are the same 
c1=unique(Hr3)
ind.location=list()
tp=1
for (i in 1:length(c1)){
  ind.loc=which(Hr3==c1[i])
  if (length(ind.loc)>1){
    ind.location[[tp]]=ind.loc
    tp=tp+1
  }

}

ind=unlist(ind.location)
MatRep_rev3=lapply(1:length(ind), function(i) Matt[[ind[i]]])

save(Hr3,ind.location,MatRep_rev3, file="checkerRevision3.RData")
#######
