
save(list=c("pheno.ori","a.indel.ori","Control.alt.counts","arity","summary.geno.extra.ori"),file="p-value.model.test.data.RData")






############# z scored based on a binomial distribution
#README
# http://en.wikipedia.org/wiki/Binomial_distribution
#  Mean in n*p
# variance np(1-p) sd =sqrt(var)
# Models use n as the coverage or reads at a given postion
## Z=(X-np)/sqrt(np(1-p)) :z-score of the minor allele frequency: (np is mean in Normal, X is measurement in tumor , the standard deviation is sqrt(np(1-p))
# and p as the frequency of the alternative allele in the control/Normal sample.

# the programs cacluales a z-score (the number of stand deviations the allele frequnecy in the tumor sample is away from the control group)
# the z-score is then converged to a p-value.
# a count for the control/tumor mox ins included which extimates the number of alternative alleles might be due to the control sample along

# cellularity is used to estimate the mix of tumor to normal tissue
# alt.reads.reference.calls is used to help estimate the alternative allele frequency in control or normal samples.
#http://www.regentsprep.org/regents/math/algtrig/ats7/blesson3.htm
#p= avg(alt reads/(normal+ alt reads))
#n= ref reads
#x= alt reads 
#abs(c)=abs(z)

p=0.5
n=100

n*p

var<- n*p*(1-p)
var

X<-80

Z<-(X-n*p)/sqrt(n*p*(1-p))
Z
#z<- (true.sample.minor.counts- true.sample.cov*normal.allele.freq) / (sqrt(true.sample.cov*normal.allele.freq*(1-normal.allele.freq)))


2*(pnorm(Z, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))

sqrt(var)
############################# rescue somatics with few calls with possion model
## cellularity<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Sequenza/NSLM.cellularity.summary.use.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)

## cellularity[1:5,]




#code.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
#working.directory<-("/media/UQCCG/UQCCG-Projects/Achal") 
working.directory<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/" # directory where   "p-value.model.test.data.RData" located  
setwd(working.directory)
source("annotate_SNPs_subroutines.r")
options("width"=250,"max.print"=1000)
core.ann<-c("chr","start","end","REF","ALT","TYPE")

#working.directory<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal" # directory where   "p-value.model.test.data.RData" located  
#setwd(working.directory)

#project.file<-"2015-03-25_SKDPNexteraRun_CSV.THP.wanted.All-maf-filtered_PAul.csv"
project.file<-"2015-03-25_SKDPNexteraRun.THP.wanted.All-maf-filtered.txt"

cellularity<-read.delim("cellularity.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)
cellularity[1:5,]

column.labels<-read.delim(project.file,header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
num.vars<-dim(column.labels)[2]
a.indel<-scan(project.file,what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
num.lines<-length(a.indel)/(num.vars)
dim(a.indel)<-c(num.vars,num.lines)
a.indel<-t(a.indel)
colnames(a.indel)<-column.labels
########################################
the.chr<-unique(a.indel[,"chr"])
print(paste("Doing Chromosome ",the.chr))

if(sum(!grepl("^chr",the.chr))>0){
  a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
}

a.indel[1:5,"chr"]
key<-build.key(a.indel,core.ann)

summary.geno.extra<-a.indel
rownames(summary.geno.extra)<-key
rownames(a.indel)<-key

#####################
#load("p-value.model.test.data.RData")

#ls() ## see what you loaded
# a.indel<-a.indel.ori
# summary.geno.extra<-summary.geno.extra.ori

#normals<-pheno.ori[pheno.ori[,"normal"],"SAMPLE"]
#normals<-normals[normals!="LPH-001-27_PD"]
#normals

Control<-c("THP-1-normal-parent")
Control
#Control.alt.counts <- alt.reads.reference.calls(a.indel,Control,AD.extension="bam.AD",threshold=1) 
Control.alt.counts<-alt.reads.Non.reference.calls(a.indel,Control,AD.extension="bam.AD",threshold=1) # 0/1 start

treatments<-c("THP-1-AS703026resist","THP-1-Selumetinibresist")
treatments


########################################
########################################
########################################
0/0 <->  0/0 P value different
0/1 <->  0/1 



snp.only<-grepl("^snp",a.indel[,"TYPE"]) ### the p.values codes uses "AD" so does ONLY work for SNPs

#has.one.geno<-as.numeric(a.indel[,"ALT.Alleles.cancer"])>0

geno.p<-genotype.p.values(a.indel[snp.only,],treatments,Control.alt.counts[snp.only, "Read.Balance"]/100,cellularity,AD.extension="bam.AD")
geno.p

#write.table(geno.p,file="/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/geno.p_0_0_nonref_0_0.txt",col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE, na="")

#p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold=2*(pnorm(abs(c(6)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold

#geno.p[geno.p>p.threshold]<-"mor.thold"

geno.p <- cbind(geno.p,Control.alt.counts[snp.only, ])
geno.p[1:10,]

snp.only.indel <- (a.indel[snp.only,])

final <- cbind (snp.only.indel[,1:92], geno.p, snp.only.indel[,93:dim(snp.only.indel)[2]])

#final<-final[!grepl("flat$",final[,"TYPE"]),]


write.table(final,file="/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/0_1_ref_0_1_with_flat.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE, na="")

found.genotype <- geno.p<=p.threshold
#geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
#geno.p[!found.genotype]<-"0/0"
#colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")













