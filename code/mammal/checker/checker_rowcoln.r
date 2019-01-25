load("ecology.RData")
source("Hfunction.r")
source("checkerBoard.r")

###########
rand1=rnorm(length(mam2))
iter=1.0e+6
#####
Hrowcol1=c()
Hrowcol1[1]=h1(mam2,rand1)
ii=c();ii1=c()
jj=c();jj1=c()
Matt=list()
Matt[[1]]=mam2
ii[1]=1;ii1[1]=1
jj[1]=2;jj1[1]=2
######
rowblock=c(1,6,29)
colblock=c(1,21,27)
#######
for (k in 2:iter){
  if (k%%2){
    TMP=CheckerBoard_revision2(Matt[[k-1]],
                               ii[k-1],ii1[k-1],jj[k-1],jj1[k-1],
                               rowBlock=rowblock,colBlock=NULL)
  }
  else 
    TMP=CheckerBoard_revision2(Matt[[k-1]],
                               ii[k-1],ii1[k-1],jj[k-1],jj1[k-1],
                               rowBlock=NULL,colBlock=colblock)

  Matt[[k]]=TMP$Matrix
  ii[k]=TMP$i
  ii1[k]=TMP$i1
  jj[k]=TMP$j
  jj1[k]=TMP$j1
  
  Hrowcol1[k]=h1(Matt[[k]],rand1)
}#end of for loop 


#####
#now check which two or three matrices are the same 
c1=unique(Hrowcol1)
ind.location=list()
tp=1
for (i in 1:length(c1)){
  ind.loc=which(Hrowcol1==c1[i])
  if (length(ind.loc)>1){
    ind.location[[tp]]=ind.loc
    tp=tp+1
  }

  if (i==length(c1)/2)
    print("half")
}

ind=unlist(ind.location)
MatRep_rowcol=lapply(1:length(ind), function(i) Matt[[ind[i]]])

save(Hrowcol1,ind.location,MatRep_rowcol, file="checkerRowCol.RData")
#######
