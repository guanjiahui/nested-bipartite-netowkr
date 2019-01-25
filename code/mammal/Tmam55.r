source("Bootstrapping (1).r")
load("ecology.RData")
library(bipartite)
##########
iter=200
Tmam55=numeric(iter)

for (l in 1:iter){
  sub11=Bootbinary(mam2[1:5,1:2])$Matrix
  sub12=Bootbinary(mam2[1:5,3:11])$Matrix  
  sub13=Bootbinary(mam2[1:5,12:15])$Matrix
  sub14=Bootbinary(mam2[1:5,16:20])$Matrix
  sub15=Bootbinary(mam2[1:5,21:26])$Matrix
  
  sub21=Bootbinary(mam2[6:8,1:2])$Matrix  
  sub22=Bootbinary(mam2[6:8,3:11])$Matrix 
  sub23=Bootbinary(mam2[6:8,12:15])$Matrix 
  sub24=Bootbinary(mam2[6:8,16:20])$Matrix
  sub25=Bootbinary(mam2[6:8,21:26])$Matrix
  
  sub31=Bootbinary(mam2[9:13,1:2])$Matrix  
  sub32=Bootbinary(mam2[9:13,3:11])$Matrix 
  sub33=Bootbinary(mam2[9:13,12:15])$Matrix 
  sub34=Bootbinary(mam2[9:13,16:20])$Matrix
  sub35=Bootbinary(mam2[9:13,21:26])$Matrix
  
  sub41=Bootbinary(mam2[14:24,1:2])$Matrix  
  sub42=Bootbinary(mam2[14:24,3:11])$Matrix  
  sub43=Bootbinary(mam2[14:24,12:15])$Matrix  
  sub44=Bootbinary(mam2[14:24,16:20])$Matrix
  sub45=Bootbinary(mam2[14:24,21:26])$Matrix
  
  sub51=Bootbinary(mam2[25:28,1:2])$Matrix 
  sub52=Bootbinary(mam2[25:28,3:11])$Matrix 
  sub53=Bootbinary(mam2[25:28,12:15])$Matrix 
  sub54=Bootbinary(mam2[25:28,16:20])$Matrix
  sub55=Bootbinary(mam2[25:28,21:26])$Matrix
  
  Binarytotal=rbind(cbind(sub11,sub12,sub13,sub14,sub15),cbind(sub21,sub22,sub23,sub24,sub25),
                    cbind(sub31,sub32,sub33,sub34,sub35),
                    cbind(sub41,sub42,sub43,sub44,sub45),cbind(sub51,sub52,sub53,sub54,sub55))
  Tmam55[l]=nestedness(Binarytotal)$temperature
}

save(Tmam55,file="Tmammal55.RData")