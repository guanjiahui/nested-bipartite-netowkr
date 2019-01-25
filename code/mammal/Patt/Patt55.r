source("exact.r")
a <- rowSums(mam2)  # row sums
b <- colSums(mam2)  # column sums
matrix_type <- 0  # 0: binary, 1: nonnegative integer

# Count the number of matrices with row sums a and column sums b
number <- count1(a,b,matrix_type)
# (or, to specify the filenames to be used:)
# number <- count(a,b,matrix_type,'input.txt','table.bin')
cat('Number of matrices =',number,'\n')
#####
number_of_samples <- 1000
x <- sample1(a,b,number_of_samples)
# (or, to specify the filenames to be used:)
# x <- sample(a,b,number_of_samples,'input.txt','table.bin','output.txt')
EPatterson=sapply(1:1000,function(i) GetBipEnergy(x[,,i]))
NODFPatterson=sapply(1:1000,function(i) nestednodf(x[,,i],order=FALSE)[3]$statistic[3])
# Display samples
for (i in 1:number_of_samples) {
  write.table(x[,,i], row.names = FALSE, col.names = FALSE)
  cat('\n')
}

#######
Count_Sample=function(Adj){
  a<- rowSums(Adj)
  b <-colSums(Adj)
  matrix_type<-0
  num1 <- count1(a,b,matrix_type)
  num1<-as.numeric(num1)
  if(num1!=1){
    x=sample1(a,b,num1)
  }
  else
    x=Adj
  return(list(num=num1,x=x))
}#Count_sample()


####################
count11=Count_Sample(mam2[1:5,1:2])$num
count12=Count_Sample(mam2[1:5,3:11])$num
count13=Count_Sample(mam2[1:5,12:15])$num
count14=Count_Sample(mam2[1:5,16:20])$num
count15=Count_Sample(mam2[1:5,21:26])$num

count21=Count_Sample(mam2[6:8,1:2])$num 
count22=Count_Sample(mam2[6:8,3:11])$num
count23=Count_Sample(mam2[6:8,12:15])$num
count24=Count_Sample(mam2[6:8,16:20])$num
count25=Count_Sample(mam2[6:8,21:26])$num

count31=Count_Sample(mam2[9:13,1:2])$num 
count32=Count_Sample(mam2[9:13,3:11])$num
count33=Count_Sample(mam2[9:13,12:15])$num
count34=Count_Sample(mam2[9:13,16:20])$num
count35=Count_Sample(mam2[9:13,21:26])$num

count41=Count_Sample(mam2[14:24,1:2])$num 
count42=Count_Sample(mam2[14:24,3:11])$num  
count43=Count_Sample(mam2[14:24,12:15])$num  
count44=Count_Sample(mam2[14:24,16:20])$num
count45=Count_Sample(mam2[14:24,21:26])$num

count51=Count_Sample(mam2[25:28,1:2])$num
count52=Count_Sample(mam2[25:28,3:11])$num 
count53=Count_Sample(mam2[25:28,12:15])$num
count54=Count_Sample(mam2[25:28,16:20])$num
count55=Count_Sample(mam2[25:28,21:26])$num
temp=1
for (i in 1:5){
  for (j in 1:5){
    C=as.numeric(eval(as.name(paste("count",i,j,sep=""))))
    if (C!=1){
      print(paste(i,j,C))
      temp=temp*C
    }
  }
}
number55=temp
number55/as.numeric(number)
############
  sample11=mam2[1:5,1:2]
  sample12=mam2[1:5,3:11]
  sample13=mam2[1:5,12:15]
  sample14=mam2[1:5,16:20]
  sample15=mam2[1:5,21:26]
  
  sample21=(mam2[6:8,1:2]) 
  sample22=(mam2[6:8,3:11])
  sample24=(mam2[6:8,16:20])
  sample25=(mam2[6:8,21:26])
  
  sample31=(mam2[9:13,1:2]) 
  sample32=(mam2[9:13,3:11])

  sample35=(mam2[9:13,21:26])
  
  sample41=(mam2[14:24,1:2]) 
  sample42=(mam2[14:24,3:11])  
  sample44=(mam2[14:24,16:20])
 
  
  sample51=(mam2[25:28,1:2])
  sample52=(mam2[25:28,3:11]) 
  sample53=(mam2[25:28,12:15])
  sample54=(mam2[25:28,16:20])

  
  sample23=Count_Sample(mam2[6:8,12:15])$x
  sample33=Count_Sample(mam2[9:13,12:15])$x
  sample34=Count_Sample(mam2[9:13,16:20])$x
  sample43=Count_Sample(mam2[14:24,12:15])$x
  sample45=Count_Sample(mam2[14:24,21:26])$x
  sample55=Count_Sample(mam2[25:28,21:26])$x
#######
Binary55=list()
temp=1
for (i1 in 1:count23){
  for (i2 in 1:count33){
    for(i3 in 1:count34){
      for (i4 in 1:count43){
        for (i5 in 1:count45){
          for (i6 in 1:count55){
            Binary55[[temp]]=rbind(cbind(sample11,sample12,sample13,sample14,sample15),
                                                      cbind(sample21,sample22,sample23[,,i1],sample24,sample25),
                                                      cbind(sample31,sample32,sample33[,,i2],sample34[,,i3],sample35),
                                                     cbind(sample41,sample42,sample43[,,i4],sample44,sample45[,,i5]),
                                                      cbind(sample51,sample52,sample53,sample54,sample55[,,i6]))
            temp=temp+1
            print(paste(i1,i2,i3,i4,i5,i6,temp))
          }
        }
      }
    }
  }
}

       
       
       
       
save.image(file="windowsEco.RData")
###############
#COUNT 2X2
B1=mam2[1:5,1:20]
B2=mam2[1:5,21:26]
B3=mam2[6:28,1:20]
B4=mam2[6:28,21:26]
 
B22count1=count1(rowSums(B1),colSums(B1),0)
B22count2=count1(rowSums(B2),colSums(B2),0)
B22count3=count1(rowSums(B3),colSums(B3),0)
B22count4=count1(rowSums(B4),colSums(B4),0)
number22=as.numeric(B22count1)*as.numeric(B22count3)*as.numeric(B22count4)
############
#Count 4x2
C1=mam2[1:5,1:2]
C2=mam2[1:5,3:11]
C3=mam2[1:5,12:20]
C4=mam2[1:5,21:26]
C5=mam2[6:28,1:2]
C6=mam2[6:28,3:11] 
C7=mam2[6:28,12:20]
C8=mam2[6:28,21:26]

B42count1=count1(rowSums(C1),colSums(C1),0)
B42count2=count1(rowSums(C2),colSums(C2),0)
B42count3=count1(rowSums(C3),colSums(C3),0)
B42count4=count1(rowSums(C4),colSums(C4),0)
B42count5=count1(rowSums(C5),colSums(C5),0)
B42count6=count1(rowSums(C6),colSums(C6),0)
B42count7=count1(rowSums(C7),colSums(C7),0)
B42count8=count1(rowSums(C8),colSums(C8),0)
number42=2*as.numeric(B42count6)*as.numeric(B42count7)*as.numeric(B42count4)
###############

number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller1=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))

number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller2=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))

number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller3=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))

number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller4=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))

number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller5=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))

number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller6=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))

##################
number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller7=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))


number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller8=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))

number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller9=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))

number_of_samples <- 99999
x <- sample1(a,b,number_of_samples)
Miller10=sapply(1:number_of_samples,function(i) h1(x[,,i],rand1))
###########
c1=unique(Miller6)
ind.location=list()
tp=1
for (i in 1:length(c1)){
  ind.loc=which(Miller6==c1[i])
  if (length(ind.loc)>1){
    ind.location[[tp]]=ind.loc
    tp=tp+1
  }

}
