source("Bootstrapping (1).r")
load("ecology.RData")
library(foreach)
library(iterators)
library(parallel)
nCores <- 8
registerDoParallel(nCores)
iter=20
Eforeach<-foreach(i = 1:iter) %dopar% {
  Bootbinary(mam2)
}

save(Eforeach,file="Eforeach.RData")