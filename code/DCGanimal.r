#total combine 4 bootstrap of animal phylogeny 
#energy plot for entry from DCG4 with animal phylogeny tree
plot(density(energy.DCGan))
#energy plot for entry from DCG4 with DCG4 tree
plot(density(energy.DCGan1))
#put them together
dat.DCGan=data.frame(dens = c(energy.DCGan,energy.DCGan1), lines = rep(c("DCG w/ an phy tree","DCG w/ DCG tree"), each=100))
densityplot(~dens,data=dat.DCGan,groups = lines,plot.points = FALSE, ref = TRUE,
            auto.key=list(lines=TRUE),xlab="Energy",main="BEEH animal (DCG as entry)")
##heatmap of the matrice
rownames(DCG4.dist)=rownames(BEEH)
colnames(DCG4.dist)=rownames(BEEH)
heatmap(DCG4.dist,margins = c(6, 8))
heatmap(DCG4.dist,Rowv=as.dendrogram(phy.an.tree),Colv=as.dendrogram(phy.an.tree),
        margins = c(6, 8)),main="BEEH DCG Entry w/ animal phylogeny ")
####################################################
#energy plot for entry from animal phylogeny w/ animal phylogeny Tree 
GetBipEnergy(BEEH.animal)
plot(density(energy.an))
plot(density(energy.an11))

####
#energy plot for entry from animal phylogeny w/ DCG4 Tree
plot(energy.anDCG)
plot(density(energy.anDCG))
plot(density(WbootT$Energy))
#####
#put them together to compare 
dat.DCGan1=data.frame(dens = c(WbootT$Energy,energy.an11), lines = rep(c("phy w/ DCG","phy w/ phy tree"), each=100))
densityplot(~dens,data=dat.DCGan1,groups = lines,plot.points = FALSE, ref = TRUE,
            auto.key=list(lines=TRUE),xlab="Energy",main="BEEH animal (anim phy as entry)")
#heapmap of the matrice
rownames(BEEH.phy)=rownames(BEEH)
colnames(BEEH.phy)=rownames(BEEH)
heatmap(BEEH.phy,Rowv=as.dendrogram(DCG.4),Colv=as.dendrogram(DCG.4),
        margins = c(6, 8))
heatmap(BEEH.phy,margins = c(6, 8))

###################################
#####
#smilar for plant 
#entry from plant phylogeny w/ DCG 1 Tree #
###########################################
plot(density(energy.plDCG1$Energy))

#####
heatmap(BEEH.phy.plant)#main="plant phylogeny")
heatmap(BEEH.phy.plant,Rowv=as.dendrogram(DCG.1),Colv=as.dendrogram(DCG.1))#,main="plant phylogeny w/ DCG tree")
####
a=rnorm(100)
PLPL.BE=BEEH.plpl1$Energy+10000*a
PLDCG.BE=BEEH.plEntryDCG1$Energy+10000*a
dat.DCGpl1=data.frame(dens = c(PLPL.BE,PLDCG.BE), lines = rep(c("plant phy w/ plant phy","plant phy w/ DCG tree"), c(100,100)))
densityplot(~dens,data=dat.DCGpl1,groups = lines,plot.points = FALSE, ref = TRUE,
            auto.key=list(lines=TRUE),xlab="Energy",main="BEEH plant (plant phy as entry)")

############
heatmap(DCG1.dist)#main="Coupling Geometry for plant")
heatmap(DCG1.dist,Rowv=as.dendrogram(phy.pl.tree),Colv=as.dendrogram(phy.pl.tree))#,
        main="Coup Geo w/ plant phylogeny tree")
##########
#entry from DCG 1 w/ DCG1 Tree 
PLPL=c(energy.plDCG3.7$Energy,energy.plDCG3.6$Energy)
DCG1DCG=PLPL[which(abs(PLPL)>abs((mean(PLPL)+0.8*sd(PLPL))))]
#entry from DCG 1 w/ plant phylogeny Tree 

DCGPL=c(energy.plDCG4.1$Energy,energy.plDCG4.6$Energy,energy.plDCG4.7$Energy)
#######
dat.DCGpl2=data.frame(dens = c(DCGPL,DCG1DCG), lines = rep(c("DCG1 w/ plant phy","DCG1 w/ DCG tree"), c(120,length(DCG1DCG))))
densityplot(~dens,data=dat.DCGpl2,groups = lines,plot.points = FALSE, ref = TRUE,
            auto.key=list(lines=TRUE),xlab="Energy",main="BEEH plant (DCG as entry)")