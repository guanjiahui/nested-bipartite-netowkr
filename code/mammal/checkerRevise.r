load("ecology.RData")
source("Hfunction.r")
source("checkerBoard.r")

###########
rand1=rnorm(length(mam2))
iter=1.0e+6
#####
Hr1=c()
Hr1[1]=h1(mam2,rand1)
ii=c();ii1=c()
jj=c();jj1=c()
Matt=list()
Matt[[1]]=mam2
ii[1]=1;ii1[1]=1
jj[1]=2;jj1[1]=2
######
rowblock=c(1,6,29)
colblock=c(1,21,27)
####
for (k in 2:iter){
  TMP=CheckerBoard_revision2(Matt[[k-1]],ii[k-1],ii1[k-1],jj[k-1],jj1[k-1])
  Matt[[k]]=TMP$Matrix
  ii[k]=TMP$i
  ii1[k]=TMP$i1
  jj[k]=TMP$j
  jj1[k]=TMP$j1

  Hr1[k]=h1(Matt[[k]],rand1)
  
  if(k==5.0e+5)
    print("HALF")
}

save(Hr1,file="checkerRevise.RData")
#######
