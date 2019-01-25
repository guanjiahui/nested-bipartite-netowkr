source("Bootstrapping (1).r")
source("EnergyOptim.r")
load("DCG.RData")
load("CoupGeo.RData")
###
iter=100
Energy.coarse1=numeric(iter)
Energy.fine1=numeric(iter)
Energy.noblock=numeric(iter)
for (l in 1:iter){
  sub1=Bootbinary(avif[1:32,1:3])$Matrix
  sub2=Bootbinary(avif[1:32,4:7])$Matrix  
  sub3=Bootbinary(avif[1:32,8:27])$Matrix
  
  sub4=Bootbinary(avif[33:45,1:3])$Matrix  
  sub5=Bootbinary(avif[33:45,4:7])$Matrix
  sub6=Bootbinary(avif[33:45,8:27])$Matrix
  Binarytotal=rbind(cbind(sub1,sub2,sub3),cbind(sub4,sub5,sub6))
  Binarytotal=rbind(Binarytotal,avif[46:56,1:27])
  Binarytotal=cbind(Binarytotal,avif[,28])
  Energy.coarse1[l]=GetBipEnergy(Binarytotal)
}

for (l in 1:iter){
  sub11=Bootbinary(avif[1:22,1:3])$Matrix
  sub12=Bootbinary(avif[23:26,1:3])$Matrix
  sub13=Bootbinary(avif[27:30,1:3])$Matrix
  sub14=Bootbinary(avif[31:32,1:3])$Matrix
  sub15=Bootbinary(avif[33:45,1:3])$Matrix
  
  sub21=Bootbinary(avif[1:22,4:7])$Matrix
  sub22=Bootbinary(avif[23:26,4:7])$Matrix
  sub23=Bootbinary(avif[27:30,4:7])$Matrix
  sub24=Bootbinary(avif[31:32,4:7])$Matrix
  sub25=Bootbinary(avif[33:45,4:7])$Matrix
  
  sub31=Bootbinary(avif[1:22,8:11])$Matrix
  sub32=Bootbinary(avif[23:26,8:11])$Matrix
  sub33=Bootbinary(avif[27:30,8:11])$Matrix
  sub34=Bootbinary(avif[31:32,8:11])$Matrix
  sub35=Bootbinary(avif[33:45,8:11])$Matrix
  
  sub41=Bootbinary(avif[1:22,12:17])$Matrix
  sub42=Bootbinary(avif[23:26,12:17])$Matrix
  sub43=Bootbinary(avif[27:30,12:17])$Matrix
  sub44=Bootbinary(avif[31:32,12:17])$Matrix
  sub45=Bootbinary(avif[33:45,12:17])$Matrix
  
  sub51=Bootbinary(avif[1:22,18:27])$Matrix
  sub52=Bootbinary(avif[23:26,18:27])$Matrix
  sub53=Bootbinary(avif[27:30,18:27])$Matrix
  sub54=Bootbinary(avif[31:32,18:27])$Matrix
  sub55=Bootbinary(avif[33:45,18:27])$Matrix
  
  temp1=rbind(sub11,sub12,sub13,sub14,sub15)
  temp2=rbind(sub21,sub22,sub23,sub24,sub25)
  temp3=rbind(sub31,sub32,sub33,sub34,sub35)
  temp4=rbind(sub41,sub42,sub43,sub44,sub45)
  temp5=rbind(sub51,sub52,sub53,sub54,sub45)
  
  Binarytotal=cbind(temp1,temp2,temp3,temp4,temp5)
  Binarytotal=rbind(Binarytotal,avif[46:56,1:27])
  Binarytotal=cbind(Binarytotal,avif[,28])
  Energy.fine1[l]=GetBipEnergy(Binarytotal)
}
for (l in 1:iter){
  Binary.no=Bootbinary(avif)$Matrix
  Energy.noblock[l]=GetBipEnergy(Binary.no)
}
#########
save(Energy.coarse1,Energy.fine1,Energy.noblock,file="avif.RData")