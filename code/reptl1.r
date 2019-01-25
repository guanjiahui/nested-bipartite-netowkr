source("Bootstrapping (1).r")
source("EnergyOptim.r")
load("DCG.RData")
load("CoupGeo.RData")
###
iter=100
Energy.overall=numeric(iter)
for (l in 1:iter){
  Binary.no=Bootbinary(reptl)$Matrix
  Energy.overall[l]=GetBipEnergy(Binary.no)
}
save(Energy.overall,file="reptl1.RData")