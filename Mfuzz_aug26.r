# This function is the core function for soft clustering. It groups genes based on the Euclidean distance
# and the c-means objective function which is a weighted square error function. Each gene is assigned
# a membership value between 0 and 1 for each cluster. Hence, genes can be assigned to different
# clusters in a gradual manner. This contrasts hard clustering where each gene can belong to a single
# cluster

library("Biobase")
library("Mfuzz")
#openVignette()
#install.packages("Mfuzz")
#install.packages("convert")
#as(object, "ExpressionSet")
zero.zero<-read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/pvalue_output_Aug21_0-0.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
rownames(zero.zero)<-zero.zero[,1]
final.zero_zero<-zero.zero
zero.one<-read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/pvalue_output_Aug21_0-1.txt",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
rownames(zero.one)<-zero.one[,1]
final.one_one<-zero.one

#jan13_2016_for amanda
all.cluster<-read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/Mfuzz/Final_Sep10/ALL_clusters_for_combined_clustering_Jan12_2016.csv",header=T,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
rownames(all.cluster)<-all.cluster[,1]
final.cluster<-all.cluster


zero.zero<-zero.zero[,255:257]
zero.one<-zero.one[,255:257]

zero.one<- all.cluster[,262:264]
##zero.zero<-zero.zero[,255:257]
##zero.zero<-cbind(id=zero.zero[,1],apply(zero.zero[,-1],2,FUN= function(x) -log(x=x,base=10))) #implementing negative log base 10
##zero.zero<-apply(zero.zero[,-1],2,FUN= function(x) -log(x=x,base=10)) #implementing negative log base 10
zero.zero<- apply(zero.zero,2,FUN= function(x) -log(x=x,base=10)) #implementing negative log base 10


##zero.one<-zero.one[,255:257]
##zero.one<-cbind(id=zero.one[,1],apply(zero.one[,-1],2,FUN= function(x) -log(x=x,base=10))) #implementing negative log base 10
##zero.one<-apply(zero.one[,-1],2,FUN= function(x) -log(x=x,base=10)) #implementing negative log base 10
zero.one<- apply(zero.one,2,FUN= function(x) -log(x=x,base=10))
##Last 30 individuals have random numbers from a Po(1)
#cc<-zero.zero[is.finite(rowSums(zero.zero)), ]
#cc<- zero.zero
cc<-zero.one
#cc <- apply(cc,1,function(x) {(x[as.numeric(x)>50 | is.na(x) | !is.finite(x)]<-50})
cc[ is.infinite(cc) | cc >50 | is.na(cc)] <- 50

# cc<-as.data.frame(zero.zero)
# cc[,2]<-as.numeric(as.character(cc[, 2] ))
# cc[,3]<-as.numeric(as.character(cc[, 3] ))
# cc[,4]<-as.numeric(as.character(cc[, 4] ))
# 
# cc <- sapply(cc,function(x) ifelse(x>50,50,x))
# cc <- sapply(cc,function(x) ifelse(x<0,0,x))
# apply(cc[,2:4],1,function(x) {x[as.numeric(x)>50 | is.na(cc) | !is.finite(cc)]<-50})
# cc[cc > 50 ] <-50
# cc[(cc[2:4] < 50),]
# tt<-zero.zero[!rowSums(zero.zero== "Inf") ,]
# pp<-tt[rowSums(tt<50),]
# rownames(zero.zero)[apply(zero.zero, 1, function(x) x[1] < x[2] & x[2] < x[3])]
# zero.zero[!rowSums((zero.zero== "Inf") | zero.zero< 50),]
##Create an expressionSet object
tmp_expr = new("ExpressionSet", exprs=cc)
mestimate(tmp_expr)

#partition coeff:Introduced by Bezdek (1981), the partition coefficient F is defined as the sum of squares of values
# of the partition matrix divided by the number of values. It is maximal if the partition is hard and
# reaches a minimum for U=1/c when every gene is equally assigned to every cluster.
# It is well-known that the partition coefficient tends to decrease monotonically with increasing n. To
# reduce this tendency we defined a normalized partition coefficient where the partition for uniform
# partitions are subtracted from the actual partition coefficients (Futschik and Kasabov,2002).

################################################ or can estimate this way
MA.expt.s.RNDM <- randomise(tmp_expr)   # MA.expt.s.RNDM <- randomise(MA.expt.s[posn.probes,])

#tmp  <- partcoef(MA.expt.s.RNDM) # This might take some time.

#tmp  <- partcoef(MA.expt.s.RNDM,crange=seq(2,6,1),mrange=seq(2.5,4.0,0.1)) #crange: range of clusters c; mrange: range of clustering parameter m.
tmp  <- partcoef(MA.expt.s.RNDM,crange=seq(2,9,1),mrange=seq(1.5,4.0,0.1)) #crange: range of clusters c; mrange: range of clustering parameter m.
#cselection(MA.expt.s.RNDM,m=3.62,crange=seq(5,40,5),repeats=5,visu=TRUE)
# c:5 c:10 c:15 c:20 c:25 c:30 c:35 c:40
# repeats:1   5    9    9   14   15   19   27   26
# repeats:2   2    9    8   15   15   16   26   23
# repeats:3   2    3   10   13   16   21   27   30
# repeats:4   5    3   11   14   19   21   29   25
# repeats:5   2    8   10   15   12   16   26   26

F<-tmp[[1]]
F.n<-tmp[[2]]
F.min<-tmp[[3]]
F>3.8*F.min
###############################################
#cl<-mfuzz(tmp_expr,c=6,m=4)
#cl2<-mfuzz(tmp_expr,c=10,m=1.5)
#cl2<-mfuzz(tmp_expr,c=15,m=2)
#cl2<-mfuzz(tmp_expr,c=8,m=1.9)
cl2<-mfuzz(tmp_expr,c=9,m=1.8) ##best
cl2<-mfuzz(tmp_expr,c=4,m=1.8) ##getting only 6 clusters
#cl2<-mfuzz(tmp_expr,c=12,m=4.1)
#mfuzz.plot(tmp_expr,cl=cl,mfrow=c(4,4),time.labels =seq(0,160,10))
#mfuzz.plot(tmp_expr,cl=cl,mfrow=c(4,4),time.labels =seq(1,50,1))
##Specify c=3 clusters
#cl = mfuzz(tmp_expr, c=3, m=1.25)
#mfuzz.plot(tmp_expr,cl=cl,mfrow=c(4,4),xlab="Samples",ylab="Pvalue(-log10)",time.labels =seq(1,50,1))
#over.lap <- overlap(cl2)
png("cluster_from_all_combined_c_9_m1.8_jan13_2016.png",
    width = 5*400,
    height = 5*400,
    res = 400,
    pointsize = 8,
    bg = "transparent")
mfuzz.plot2(tmp_expr,cl=cl2,mfrow=c(3,3),xlab="Samples",ylab="Pvalue(-log10)", x11=FALSE)
dev.off()

png("cluster_from_all_combined_c_4_m1.8_jan13_2016.png",
    width = 5*400,
    height = 5*400,
    res = 400,
    pointsize = 8,
    bg = "transparent")
mfuzz.plot2(tmp_expr,cl=cl2,mfrow=c(3,3),xlab="Samples",ylab="Pvalue(-log10)", x11=FALSE)
dev.off()


# $cluster: a vector of integers containing the indices of the clusters where the data points
# are assigned to for the closest hard clustering, as obtained by assigning points to
# the (first) class with maximal membership
membership <- cl2$cluster 
# a matrix with the membership values of the data points to the clusters.
#A soft cluster is considered as empty, if none of the genes has a corresponding membership value larger than 0.5.
confidence<-cl2$membership
confidence<-apply(confidence,1,max)
#the number of data points in each cluster of the closest hard clustering
cl2$size
cc[1:5,]
zero.zero[1:5,]
length(confidence)

dim(tmp_expr)
gg<-cbind(final.one_one,confidence,membership)
gg<-cbind(final.cluster,confidence,membership)
#gg<-cbind(final.one_one,confidence,membership)
confidence[1:5]




write.table(gg,file="member_confidence_cl2_Sep3_0-1.txt",col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)

write.table(gg,file="member_confidence_cluster_from_all_combined_c_9_m1.8_jan13_2016_cl2.txt",col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE)


#fig.prefix
#savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")


#############################################################
#############################################################

pDataFile<-ExpressionSet(assayData=zero.zero)


minimalSet <- new("ExpressionSet", exprs = exprs)

   http://itb.biologie.hu-berlin.de/futschik/software/R/Mfuzz
   
   data(yeast)
data(zero.zero)
   yeast.r<-filter.NA(yeast,thres=0.25)

   yeast.f<-fill.NA(yeast.r,mode="mean")             # mode =knn,wknn
   
   yeast.s<-standardise(yeast.f)
   
   class(yeast.s)
   class(zero.zero)
   cl<-mfuzz(yeast.s,c=16,m=1.25)
mfuzz.plot(yeast.s,cl=cl,mfrow=c(4,4),time.labels=seq(0,160,10))


   # cluster stability:
cl2 <- mfuzz(yeast.s, c = 16, m = 1.35)   #M->1 is hard clustering  0 or 1 membership .
                        # expression vectots with large noise have  low membership so large M gives measure of noise
                         # m controls sensitivity to noise
 partcoef                 # empty cluster all membershipd < 0.5

  #### parameter selection
yeastFR <- randomise(yeast.s)
cl <- mfuzz(yeastFR,c=20,m=1.1)
mfuzz.plot(yeastFR,cl=cl,mfrow=c(4,5)) # shows cluster structures (non-uniform partition)

cl$centers

 tmp  <- partcoef(yeastFR) # This might take some time.
 F <- tmp[[1]];F.n <- tmp[[2]];F.min <- tmp[[3]]

 # Which clustering parameters result in a uniform partition?  
 F > 1.01 * F.min

cl <- mfuzz(yeastFR,c=20,m=1.25) # produces uniform partion      # munifrom color is uniform coloring

mfuzz.plot(yeastFR,cl=cl,mfrow=c(4,5))
# uniform coloring of temporal profiles indicates uniform partition
}


 tmp  <- cselection(yeast.s,m=1.25,crange=seq(5,40,5),repeats=5,visu=TRUE)

 

                         

par(mfrow = c(1,1))
mfuzz.plot(yeast.s, cl = cl2, mfrow = c(4, 4), time.labels = seq(0,
+ 160, 10))



O <- overlap(cl)
Ptmp <- overlap.plot(cl, over = O, thres = 0.05)


   > cl3 <- mfuzz(yeast.s, c = 10, m = 1.25)
> mfuzz.plot(yeast.s, cl = cl3, mfrow = c(3, 4))
> O3 <- overlap(cl3)
> overlap.plot(cl3, over = O3, P = Ptmp, thres = 0.05)


   > cl4 <- mfuzz(yeast.s, c = 25, m = 1.25)
> mfuzz.plot(yeast.s, cl = cl4, mfrow = c(5, 5))
> O4 <- overlap(cl4)
> overlap.plot(cl4, over = O4, P = Ptmp, thres = 0.05)



###################################################


cl2<-mfuzz(MA.expt.s.RNDM,c=8,m=3.0)   # cl2<-mfuzz(MA.expt.s[posn.probes,],c=10,m=1.26)
mfuzz.plot(MA.expt.s.RNDM,cl=cl2,mfrow=c(3,4))

################################

cl.1<-mfuzz(MA.expt.s1,c=8,m=2.5,iter.max=5000)
#cl2<-mfuzz(MA.expt.s,c=6,m=2.05,iter.max=500)  # For carcin # cl2<-mfuzz(MA.expt.s[posn.probes,],c=10,m=1.26)

#cl2<-mfuzz(MA.expt.s,c=6,m=2.2)


mfuzz.plot(MA.expt.s,cl=cl2,mfrow=c(3,3),time.labels=c(3,24,27,32,48),min.mem=0.5)

membership<-cl2$cluster # cloest hard cluster
confidence<-cl2$membership
confidence<-apply(confidence,1,max)
cl2$size



###########################################################################################
library(Mfuzz)

tps = 6;cases = 90
d = rnorm(tps*cases, 1)  ##Poisson distribution with mean 1
m = matrix(d, ncol=tps, nrow=cases)

##First 30 individuals have increasing trends
m[1:30,] = t(apply(m[1:30,], 1, cumsum))

##Next 30 have decreasing trends
##A bit hacky, sorry
m[31:60,] = t(apply(t(apply(m[31:60,], 1, cumsum)), 1, rev))

##Last 30 individuals have random numbers from a Po(1)

##Create an expressionSet object
tmp_expr = new('ExpressionSet', exprs=m)

##Specify c=3 clusters
cl = mfuzz(tmp_expr, c=6, m=1.25)
mfuzz.plot(tmp_expr,cl=cl, mfrow=c(2, 2))


###########################################################################################





fig.prefix
savePlot(filename=paste(fig.prefix,".png",sep=''),type="png")
savePlot(filename=paste(fig.prefix,".jpeg",sep=''),type="jpeg")
savePlot(filename=paste(fig.prefix,".tiff",sep=''),type="tiff")
#savePlot(filename=paste(fig.prefix,"aa.bmp",sep=''),type="bmp")
dev.print(svg,paste(fig.prefix,".svg",sep=''))
