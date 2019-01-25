load("ecology.RData")
source("Hfunction.r")
source("checkerBoard.r")

###########
rand1=rnorm(length(mam2))
iter=1.0e+6
#####
Hr2=c()
Hr2[1]=h1(mam2,rand1)
ii=c();ii1=c()
jj=c();jj1=c()
Matt=list()
Matt[[1]]=mam2
ii[1]=1;ii1[1]=1
jj[1]=2;jj1[1]=2
######
for (k in 2:iter){
   TMP=CheckerBoard_revision3(Matt[[k-1]],
                               ii[k-1],ii1[k-1],jj[k-1],jj1[k-1])

  Matt[[k]]=TMP$Matrix
  ii[k]=TMP$i
  ii1[k]=TMP$i1
  jj[k]=TMP$j
  jj1[k]=TMP$j1
  
  Hr2[k]=h1(Matt[[k]],rand1)
}#end of for loop 

###########
#now check which two or three matrices are the same 
c1=unique(Hr2)
ind.location=list()
tp=1
for (i in 1:length(c1)){
  ind.loc=which(Hr2==c1[i])
  if (length(ind.loc)>1){
    ind.location[[tp]]=ind.loc
    tp=tp+1
  }

  if (i==length(c1)/2)
    print("half")
}

ind=unlist(ind.location)
MatRep_rev2=lapply(1:length(ind), function(i) Matt[[ind[i]]])

save(Hr2,ind.location,MatRep_rev2, file="checkerRevision2.RData")
#######


############
#combined algorithm: Revised checker-board +Miller's 
Energy_rev2=c()
Energy_rev2[1]=GetBipEnergy(mam2)
ii=c();ii1=c()
jj=c();jj1=c()
Matt=list()
Matt[[1]]=mam2
ii[1]=1;ii1[1]=1
jj[1]=2;jj1[1]=2
######
iter=10000
for (k in 2:iter){
  TMP=CheckerBoard_revision3(Matt[[k-1]],
                             ii[k-1],ii1[k-1],jj[k-1],jj1[k-1])
  
  Matt[[k]]=TMP$Matrix
  ii[k]=TMP$i
  ii1[k]=TMP$i1
  jj[k]=TMP$j
  jj1[k]=TMP$j1
  
  Energy_rev2[k]=GetBipEnergy(Matt[[k]])
}#end of for loop 
NODF_rev2=sapply(1:10000,function(i) nestednodf(Matt[[i]],order=FALSE)[3]$statistic[3])
T_rev2=sapply(1:10000,function(i) nestedtemp(Matt[[i]])$statistic)
Matt_temp=lapply(1:10000, function(i) t(Matt[[i]])[rev(1:ncol(Matt[[i]])),])
N_rev2=sapply(1:10000,function(i) Nplus(Matt_temp[[i]]))
############

#next try use the Miller's algorithmn
number_of_samples <- 10000
x <- sample1(a,b,number_of_samples)
EPatterson=sapply(1:10000,function(i) GetBipEnergy(x[,,i]))
NODFPatterson=sapply(1:10000,function(i) nestednodf(x[,,i],order=FALSE)[3]$statistic[3])
TPatterson=sapply(1:10000,function(i) nestedtemp(x[,,i])$statistic)
Patter_temp=lapply(1:10000, function(i) t(x[,,i])[rev(1:ncol(x[,,i])),])
N_Patterson=sapply(1:10000,function(i) Nplus(Patter_temp[[i]]))
###########

# Next do the checker-Board Swapping again
Energy2_rev2=c()
Energy2_rev2[1]=GetBipEnergy(Matt[[1000]])
ii=c();ii1=c()
jj=c();jj1=c()
Matt2=list()
Matt2[[1]]=Matt[[1000]]
ii[1]=1;ii1[1]=1
jj[1]=2;jj1[1]=2
##
iter=10000
for (k in 2:iter){
  TMP=CheckerBoard_revision3(Matt2[[k-1]],
                             ii[k-1],ii1[k-1],jj[k-1],jj1[k-1])
  
  Matt2[[k]]=TMP$Matrix
  ii[k]=TMP$i
  ii1[k]=TMP$i1
  jj[k]=TMP$j
  jj1[k]=TMP$j1
  
  Energy2_rev2[k]=GetBipEnergy(Matt2[[k]])
}#end of for loop 
NODF2_rev2=sapply(1:10000,function(i) nestednodf(Matt2[[i]],order=FALSE)[3]$statistic[3])
T2_rev2=sapply(1:10000,function(i) nestedtemp(Matt2[[i]])$statistic)
Matt_temp2=lapply(1:10000, function(i) t(Matt2[[i]])[rev(1:ncol(Matt2[[i]])),])
N2_rev2=sapply(1:10000,function(i) Nplus(Matt_temp2[[i]]))
############

# check the distribution of different index
#Energy
dat.E<- data.frame(dens = c(Energy2_rev2,EPatterson,ENature), 
                   lines = rep(c("Checker-board","Miller's","Nature Communication"),
                               c(length(Energy2_rev2),length(EPatterson),length(ENature))))
densityplot(~dens,data=dat.E,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="Energy Comparison")

#NODF
dat.NODF<- data.frame(dens = c(NODF2_rev2,NODFPatterson,NODFNature), 
                   lines = rep(c("Checker-board","Miller's","Nature Communication"),
                               c(length(NODF2_rev2),length(NODFPatterson),length(NODFNature))))
densityplot(~dens,data=dat.NODF,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="NODF Comparison")

#Temperature
dat.T<- data.frame(dens = c(T2_rev2,TPatterson,TNature), 
                      lines = rep(c("Checker-board","Miller's","Nature Communication"),
                                  c(length(T2_rev2),length(TPatterson),length(TNature))))
densityplot(~dens,data=dat.T,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="Temperature Comparison")

#Nplus
dat.N<- data.frame(dens = c(N2_rev2,N_Patterson,NplusNature), 
                   lines = rep(c("Checker-board","Miller's","Nature Communication"),
                               c(length(N2_rev2),length(N_Patterson),length(NplusNature))))
densityplot(~dens,data=dat.N,groups = lines,plot.points = FALSE, 
            ref = TRUE, xlab="Energy",
            auto.key=list(lines=TRUE),main="Nplus Comparison")
####################
plot(c(Energy2_rev2,EPatterson,ENature),type="l",main="Energy")
plot(c(NODF2_rev2,NODFPatterson,NODFNature),type="l",main="NODF")
plot( c(T2_rev2,TPatterson,TNature),type="l",main="Temperature")
plot(c(N2_rev2,N_Patterson,NplusNature),type="l",main="Nplus")
