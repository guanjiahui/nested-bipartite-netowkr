source("Bootstrapping (1).r")
load("ecology.RData")
library(bipartite)
##########
iter=200

Matrix42=list()
for (l in 1:iter){
  sub1=mam2[1:5,1:2]
  sub2=Bootbinary(mam2[1:5,3:11])$Matrix 
  sub3=Bootbinary(mam2[1:5,12:20])$Matrix
  sub4=Bootbinary(mam2[1:5,21:26])$Matrix
  
  sub5=mam2[6:28,1:2]
  sub6=Bootbinary(mam2[6:28,3:11])$Matrix 
  sub7=Bootbinary(mam2[6:28,12:20])$Matrix
  sub8=Bootbinary(mam2[6:28,21:26])$Matrix
  
  Binarytotal=rbind(cbind(sub1,sub2,sub3,sub4),cbind(sub5,sub6,sub7,sub8))
  Matrix42[[l]]=Binarytotal
}

save(Matrix42,file="Matrix42.RData")

#######################################
iter=5000

sub1=mam2[1:5,1:2]
sub2=Miller_regular(mam2[1:5,3:11],iter)$Matrix 
sub3=Miller_regular(mam2[1:5,12:20],iter)$Matrix
sub4=Miller_regular(mam2[1:5,21:26],iter)$Matrix

sub5=mam2[6:28,1:2]
sub6=Miller_regular(mam2[6:28,3:11],iter)$Matrix 
sub7=Miller_regular(mam2[6:28,12:20],iter)$Matrix
sub8=Miller_regular(mam2[6:28,21:26],iter)$Matrix


#E.Miller.42=c()
#TMiller.42=c()
#NODF.Miller.42=c()
#Nplus.Miller.42=c()
#DM.42=c()
#DM2.42=c()

NCG42=c()


#Checker_Diamond42=c()
#Cscore_Stone42=c()
#Vratio_Roboson42=c()
#Combo_Pielou42=c()

for (l in 1:iter){
  Binarytotal=rbind(cbind(sub1,sub2[,,l],sub3[,,l],sub4[,,l]),
                    cbind(sub5,sub6[,,l],sub7[,,l],sub8[,,l]))
  #E.Miller.42[l]=GetBipEnergy(Binarytotal)
  #NODF.Miller.42[l]=nestednodf(Binarytotal,order=FALSE)[3]$statistic[3]
  #TMiller.42[l]=nestedtemp(Binarytotal)$statistic
 # Nplus.Miller.42[l]=Nplus( t(Binarytotal)[rev(1:ncol(Binarytotal)),])
  #DM.42[l]=DM(t(Binarytotal),c,r,5,5)
  #DM2.42[l]=DM4(t(Binarytotal),c1,r1,5,5,t(mam3))
  #GG42[l]=GG(t(Binarytotal),c1,r1,rowchoose,5,5)
  NCG42[l]=N_CG(t(Binarytotal),c1,r1,5,5,t(mam3))
  
  
 # Checker_Diamond42[l]=checker(Binarytotal)
#  Cscore_Stone42[l]= c_score(Binarytotal)
#  Vratio_Roboson42[l]=v_ratio(Binarytotal)
#  Combo_Pielou42[l]=COMBO(Binarytotal)
}