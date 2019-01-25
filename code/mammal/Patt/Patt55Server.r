lowsource("windowsEco.RData")
source("ecology.RData")
Binary55=list()
temp=1
for (i1 in 1:count23){
  for (i2 in 1:count33){
    for(i3 in 1:count34){

    }#for i3
  }#for i2
}#for i1
Tot=prod(count23,count33,count34,count43,count45,count55)
EPatt55=sapply(1:Tot,function(i) GetBipEnergy(Binary55[[i]]))

save(Binary55,EPatt55,file="Patt55.RData")


temp=1
temp2=1
LOW_matrix2=list()
for (i4 in 1:count43){
  for (i5 in 1:count45){
    for (i6 in 1:count55){
      Binary55[[temp]]=rbind(cbind(sample11,sample12,sample13,sample14,sample15),
                             cbind(sample21,sample22,sample23[,,4],sample24,sample25),
                             cbind(sample31,sample32,sample33[,,4],sample34[,,2],sample35),
                             cbind(sample41,sample42,sample43[,,i4],sample44,sample45[,,i5]),
                             cbind(sample51,sample52,sample53,sample54,sample55[,,i6]))
      E_check=GetBipEnergy(Binary55[[temp]])
      if (E_check<=E_lowest){
        E_lowest=E_check
        LOW_matrix2[[temp2]]=Binary55[[temp]]
        temp2=temp2+1
      }
      temp=temp+1
    }#for i6
  }#for i5
  print(i4)
}#for i4

E_find_low=sapply(82620:(temp-1), function(i) GetBipEnergy(Binary55[[i]]))
cc=which(E_find_low==-2204)
low_matrix1=list()
for (i in 1:length(cc)){
  low_matrix1[[i]]=Binary55[[cc[i]]]
}
E_lowest=-2204
#####
indunique=c(1,7,9,11,13,19,21,23,26,31,33,36,38,43,45,47)

nodf_low2=sapply(1:16,function(i) nestednodf(LOW_matrix2[[indunique[i]]],order=FALSE)[3]$statistic[3])
T_low2=sapply(1:16,function(i) nestedtemp(LOW_matrix2[[indunique[i]]])$statistic)
N_low2=sapply(1:16,function(i) Nplus(LOW_matrix2[[indunique[i]]]))
#############
c1=unique(H222)
ind.location=list()
tp=1
for (i in 1:length(c1)){
   ind.loc=which(H222==c1[i])
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
