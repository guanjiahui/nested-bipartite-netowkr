source("Bootstrapping (1).r")
load("ecology.RData")

##########
col.order2=c(2,1,4,3,5,6,7,8,9,16,17,11,10,12,13,14,19,18,15,20,22,21,23,24,25,26)
mam3=t(mam)[row.order,col.order2]
A=mam3[6:28,12:16]
B=mam3[6:28,17:20]
AB=mam3[6:28,12:20]
########
iter=200
E.AandB2=numeric(iter)
EAB2=numeric(iter)


for (l in 1:iter){
  sub1=Bootbinary(A)$Matrix
  sub2=Bootbinary(B)$Matrix  
  
  Binarytotal=cbind(sub1,sub2)
  E.AandB2[l]=GetBipEnergy(Binarytotal)
  #####
  EAB2=Bootbinary(AB)$Energy
}

save(E.AandB2,EAB2,file="AB2.RData")