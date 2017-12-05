library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(XLConnect)
library("openxlsx")
library("data.table")
library("dplyr")
#library("plyr")

code.dir <- "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
core.ann <- c( "chr","start","end","ALT","REF","TYPE")
core.ann.gr <- c( "chr","start","end","strand")
#mydf  <-  read.delim("mydf_enrichment_test_file.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="")
mydf  <-  read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/2015-03-25_SKDPNexteraRun.THP.wanted.All-maf-filtered.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,strip.white=TRUE,na.strings="",quote="")
finaldf<-mydf  # to cbind the final coverage output
a.indel.mydf<-as.matrix(finaldf)

#Triming leading and trailing whitespace
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

a.indel.mydf <-trim(a.indel.mydf)
a.indel.mydf[1:5,1:5]
###########################
chr.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-25_SKDPNexteraRun/Analysis/"
all.chr.files<-list.files(chr.dir)
all.chr.files<-all.chr.files[grepl("ALL.ALL_GENOTYPES_analysis.txt",all.chr.files)]
output.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/wanted_snps/Analysis/"

#this loop below will find all the matching snps present in mydf and in all the chomosomes 
i<-1
for(i in 1:length(all.chr.files)){
setwd("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-25_SKDPNexteraRun/Analysis")
#a.indel from geno dump
#a.indel<-read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-25_SKDPNexteraRun/Analysis/2015-03-25_SKDPNexteraRun.chr2.ALL.ALL_GENOTYPES_analysis.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,strip.white=TRUE,na.strings="",quote="")
a.indel<- read.delim(all.chr.files[i],header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,strip.white=TRUE,na.strings="",quote="")

#project.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-25_SKDPNexteraRun/Analysis/2015-03-25_SKDPNexteraRun.chr2.ALL.ALL_GENOTYPES_analysis.txt"
############################################################
#project.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/2013-10-27_AML_with_AOGCControl_NoFailedLane.chr13.ALL.ALL_GENOTYPES_WITH_SYNON_analysis.txt"
##################################################
#project.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/2013-10-27_AML_with_AOGCControl_NoFailedLane.chr13_predictions.txt"
#project.file<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2013/2013-10-27_AML_with_AOGCControl_NoFailedLane/Analysis/2013-10-27_AML_with_AOGCControl_NoFailedLane.chr13.TGCM-AML.analysis.txt"
##################################################
#setwd(analysis.dir)
###################################
a.indel[1:5,1:5]
key.indel<-build.key(a.indel,core.ann)
key.mydf<-build.key(mydf,core.ann)
wanted.pos<-match(key.indel,key.mydf)
missing<-is.na(wanted.pos)
wanted.snps<-a.indel[!missing,]
setwd(output.dir)
write.table(wanted.snps,file=all.chr.files[i],col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE, na="")
rm(a.indel)
}

######################################
#this wil merge all the snps in 24 chromosomes

UQCCG.data<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/" #directoy equivalent to UQCCG.data 
project<-"wanted_snps"
project.name <- project
#annotate.dir <-paste(UQCCG.data,project,"Annotate",sep="/")
analysis.dir<-paste(UQCCG.data,project,"Analysis",sep="/")

output.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/wanted_snps/merged/"


#project.extension<-"analysis.txt" # use this one for compete listing
#project.extension<-"variants.analysis.txt" ; fam<-c(".large.")
#project.extension<-"analysis-maf-filtered.txt"
project.extension<-"analysis.txt"
#".analysis-maf-filtered.txt"   ".analysis-maf-filtered.txt.geno.all.txt"  ## just the exterion not fam.extension! ".analysis-maf-filtered.txt.geno.all.txt" #
fam<-c(".ALL.ALL_GENOTYPES_") #  ALL or  all ""-one project (the prefix of the summary files to collect
#the.sample.sheet<-"/mnt/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-03-19_MND_MODY_LKAF_Nimbelgen/BAM/140206_SN775_0138_AC3751ACXX-SampleSheet_MODY_SJ.csv"
the.sample.sheet<-"no_file" # Use no_file to just use all in the genotyping run "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2014-06-24_AML_RSGB_AOGC_withHaplotypeCaller/BAM/TGCM-AML_RSGB_PILOT_SampleSheet_controls.csv"

remove.from.controls<-c() # expand.labels.to.samples(remove.from.controls,control.samples)
remove.from.all.samples<-c() #expand.labels.to.samples(remove.from.all.samples,all.samples)

# 300 samples is too many and exceeds memory limit
# remove.from.controls<-c("Q","AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q") # expand.labels.to.samples(remove.from.controls,control.samples)
# remove.from.all.samples<-c("AD","DP","0.005(MAF).Q","0.001(MAF).Q","0(MAF).Q","Q.MAF:0.001","Q.MAF:0.005","Q.MAF:0") #expand.labels.to.samples(remove.from.all.samples,all.samples)

remove.cols<-c(     ) ### may cause error reload after read in a.indels
dont.build.summary<-TRUE # FALSE then get a combined result
accumulate<-TRUE
GQ.thresh<-0 ### this will only write FILTER=PASS
# # ###############################################


############################################# BEGIN #######################################################
options(stringsASFactors=FALSE)
options(width=250,max.print=4000)
genome.build<-"hg19"
if(!exists("dont.build.summary")){dont.build.summary<-FALSE} ## if TRUE does not make summary indel file
if(!exists("set.unlabelled.to.control")){set.unlabelled.to.control<-TRUE} # set samples not listed in the samplesheet to controls 
######################################################### Set up the basics for each run #####
########################### Set up for annovar #################################     
maf.threshold<-0.0 
if(!exists("GQ.thresh")){GQ.thresh<-0}
anno.DB.location.core<-"/mnt/UQCCG/Software/annovar/humandb"
anno.DB.location<-paste(anno.DB.location.core,genome.build,sep="/")

######################### Info the run
generic.filter.DB<-c("hg19_sequnome.txt")    # returns 2  extra columns (DB,score)
names(generic.filter.DB)<-c("AML-sequnome")
update.annovar.annotations<-TRUE
core.ann<-c("chr","start","end","REF","ALT","TYPE")
##################################
code.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts/"
setwd(code.dir)
source("ucsc.table.names.r")   # load in the UCSC tables these use the db file names and not their lable-names 
source("annotate_SNPs_subroutines.r")
source("ucsc.table.names.processor.r")





#### test fam list
files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
toString( unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   ))) )
toString( files)
####


#############################################################################################################

######################################### Predefined variables required
##################################################################################


#### assume has format project.chr.fam.extension or chr.project.fam.extension
setwd(analysis.dir)

files<-dir(analysis.dir)
the.extension<-paste(project.extension,"$",sep="")
files<-files[grepl(the.extension ,files)]
if(fam=="ALL" | fam=="All" | fam=="all" ){
  fam<-unique(unlist( mapply(function(x){x[length(x)]}, strsplit(gsub(the.extension,"",files),split=".",fixed=TRUE)   )))
}

fam

ifam<-1
for(ifam in 1:length(fam)){
  
  the.extension<-paste(fam[ifam],project.extension,"$",sep="")
  project.files<-files[grepl(the.extension ,files)]
  print(sort(paste("Doing: ",project.files,sep=""))) # project.files<-project.files[1:22]
  
  indels<-{}
  the.col<-{}
  
  
  data.summary<-{} # used if doing genotype counts
  project.files
  # ichr<-1
  
  
  
  
  
  project.files<-project.files[!grepl("chrALL",project.files)] ## remove old combined runs 
  
    
  if(length(project.files)!=24){
    print("########################################### WARNING #################################")
    print("less that 24 chromosomes detected")
    print(fam[ifam])
    print("########################################### WARNING #################################") 
  }
  project.files  #  project.files <- project.files[-23]
  ichr<-19
  for(ichr in 1:length(project.files)){
    print(project.files[ichr])
    
    
    setwd(analysis.dir)
    ################## fast read ###########
    column.labels<-read.delim(project.files[ichr],header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
    num.vars<-dim(column.labels)[2]
    a.indel<-scan(project.files[ichr],what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="") # quote="\""
    num.lines<-length(a.indel)/(num.vars)
    dim(a.indel)<-c(num.vars,num.lines)
    a.indel<-t(a.indel)
    colnames(a.indel)<-column.labels
    ########################################
    #a.indel<-read.delim(project.files[ichr],header=T,skip=0,fill=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
    all.possible.samples<-gsub(".GT$","",colnames(a.indel)[grep(".GT$",colnames(a.indel))],perl=TRUE)
   
    print(dim(a.indel))
    a.indel[1:2,1:50]
    
    #    table(a.indel[,"FILTER"])
    
    
    
    ######################### subset samples #################
    ########################################### samlpe subsetting using sample sheet
    
    ### no_file just use all samples - all assumed as affected
    if(the.sample.sheet=="no.file" | the.sample.sheet=="no_file"){
      all.samples<-all.possible.samples
      sample.sheet.full<-cbind(fam[ifam],fam[ifam], all.samples,".",".","-9",2)
      colnames(sample.sheet.full)<-c("SampleProject","FamilyCode","ParticipantCode","PaternalID","MaternalID","Sex","AffectionStatus")
      control.samples<-{} # assume no controls  control.samples<-c()
      
      
      
    }else{
      sample.sheet.full<-read.delim(the.sample.sheet,header=T,sep=",",fill=TRUE,stringsAsFactors=FALSE)
      sample.sheet.full[1:5,]
      dim(sample.sheet.full)
      colnames(sample.sheet.full)
      tapply(sample.sheet.full[,"SampleProject"],sample.sheet.full[,"SampleProject"],length)
      
      all.samples<-unique(sample.sheet.full[,"ParticipantCode"])
      control.samples<-unique(sample.sheet.full[sample.sheet.full[,"SampleProject"]=="Control" | sample.sheet.full[,"SampleProject"]=="control"  | sample.sheet.full[,"SampleProject"]=="CONTROL" , "ParticipantCode"])
    }
    
    sample.sheet.full[1:5,]
    
    #######################identify cols to remove
    
    
    
    if(!is.null(remove.from.controls) & !is.null(control.samples)){
      to.remove.control<-expand.labels.to.samples.complex(remove.from.controls,control.samples,paste.after=FALSE,seperator=".")
      to.remove.control.2<-expand.labels.to.samples.complex(remove.from.controls,control.samples,paste.after=FALSE,seperator=":") ## get "S02-F21-P01:0.005(MAF).Q"
      to.remove.control<-unique(c(to.remove.control,to.remove.control.2))
    }else{
      to.remove.control<-{}
    }
    
    if(!is.null(remove.from.all.samples)){
      to.remove.all<-expand.labels.to.samples(remove.from.all.samples,all.samples)
    }else{
      to.remove.all<-{}
    }
    
    to.remove.samples<-unique(c(to.remove.control,to.remove.all))
    remove.cols<-unique(c(remove.cols,to.remove.samples))
    
    
    if(!is.null(remove.cols)){
      a.indel<-a.indel[,colnames(a.indel)[!(colnames(a.indel) %in% remove.cols)]]
    }
   
    
    if(is.null(dim(indels))){
      indels<-a.indel
      the.col<-colnames(a.indel)
    }else{
      
      if(sum(!(the.col %in% colnames(a.indel)))>0  | sum(!(colnames(a.indel) %in% the.col))>0){
        print("error colnames don't match")
        print(the.col[!(the.col %in% colnames(a.indel))])
        print(colnames(a.indel)[!(colnames(a.indel) %in% the.col)])
        next
      } # columns don't match
      
      if(sum(!(colnames(a.indel) %in% the.col))>0 ){
        a.indel<-a.indel[,the.col]     } # reduce a.indels }
      
      indels<-rbind(indels,a.indel[,the.col])
    } ## is null so  not first
    #print(project.files[ichr])
    #print(dim(a.indel))
    
    # HERE MHAIRI 
        
    ############################ GQ.thresh<-20
    if(accumulate){
      
      setwd(output.dir)
      #   best.quality<-as.numeric(a.indel[,"GQ_MEAN"])>=GQ.thresh & a.indel[,"FILTER"]=="PASS"
      # best.quality<- a.indel[, "FILTER"]=="PASS" #  &  a.indel[,"wanted.muts"]=="TRUE" #  & a.indel[,"MAF.lt:0.005"]=="TRUE"
      best.quality<-a.indel[,"MAF.lt:0.025"]=="TRUE" | a.indel[,"Gene.Names"] %in% c("DDX41","TET2", "GATA2", "ASXL1", "NOTCH1", "IDH1", "JAK2","WT1","MLL","KRAS","FLT3","IDH2","IDH1","TP53","KIT","NPM1","JAK2","DNMT3A","TET2","RUNX1","NRAS","CEBPA","PTPN11","U2AF1","SMC1A","SMC3","PHF6","STAG2","RAD21","FAM5C","EZH2","HNRNPK","FANCA","FANCB","FANCC","FANCD1","FANCD2","FANCE","FANCF","FANCG","FANCI","BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4","ERCC4","APITD1","STRA13","C1orf86","C19orf40","C17orf70","SLX1","MUS81","ERCC1","FAN1","EME1","EME2","MRE11A","NBN1","RAD50","FAND1","BRCA1","BARD1","RAD51","RAD51B","RAD51D","XRCC2","XRCC3","RMI1","RMI2","BLM","TOP3A","RPA1","RPA2","RPA3","ATM","ATR","ATRIP","CHECK1","RAD9A","RAD17","CS","DLAT","DLD","DLST","FH","IDH1","IDH2","IDH3A","IDH3B","IDH3G","MDH1","MDH2","ACLY","ACO1","OGDH","ACO2","PC","PCK1","PCK2","PDHA1","PDHA2","PDHB","OGDHL","SDHA","SDHB","SDHC","SDHD","SUCLG2","SUCLG1","SUCLA2")
      
      best.quality2<-a.indel[,"Gene.Names"] %in%  c("DDX41","TET2", "GATA2", "ASXL1", "NOTCH1", "IDH1", "JAK2","WT1","MLL","KRAS","FLT3","IDH2","IDH1","TP53","KIT","NPM1","JAK2","DNMT3A","TET2","RUNX1","NRAS","CEBPA","PTPN11","U2AF1","SMC1A","SMC3","PHF6","STAG2","RAD21","FAM5C","EZH2","HNRNPK","FANCA","FANCB","FANCC","FANCD1","FANCD2","FANCE","FANCF","FANCG","FANCI","BRIP1","FANCL","FANCM","PALB2","RAD51C","SLX4","ERCC4","APITD1","STRA13","C1orf86","C19orf40","C17orf70","SLX1","MUS81","ERCC1","FAN1","EME1","EME2","MRE11A","NBN1","RAD50","FAND1","BRCA1","BARD1","RAD51","RAD51B","RAD51D","XRCC2","XRCC3","RMI1","RMI2","BLM","TOP3A","RPA1","RPA2","RPA3","ATM","ATR","ATRIP","CHECK1","RAD9A","RAD17","CS","DLAT","DLD","DLST","FH","IDH1","IDH2","IDH3A","IDH3B","IDH3G","MDH1","MDH2","ACLY","ACO1","OGDH","ACO2","PC","PCK1","PCK2","PDHA1","PDHA2","PDHB","OGDHL","SDHA","SDHB","SDHC","SDHD","SUCLG2","SUCLG1","SUCLA2")
      
      
      #      best.quality <-best.quality # & !is.na(best.quality)
      sum(best.quality)
      
      
      if(ichr==1){
        max.maf<-max(as.numeric(gsub("^MAF.lt:","",colnames(a.indel)[grepl("^MAF.lt:",colnames(a.indel))])))
        
        max.maf.col<-paste("MAF.lt:",max.maf,sep="")
      }
      
      #     an.indel<-grepl("^indel",a.indel[,"TYPE"])
      #sum(an.indel)
      ##      dim(a.indel[an.indel,])
      
      if(ichr==1){
        write.table(a.indel,file=paste(project.name,".BEST.chrALL.ACC_ACHEL",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
        # write.table(a.indel[best.quality,],file=paste(project.name,".BEST.chrALL.ACC_ACHEL_0.025",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
        #write.table(a.indel[best.quality2,],file=paste(project.name,".BEST.chrALL.ACC_ACHEL_SUBSET",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
        ## write.table(a.indel[an.indel,],file=paste(project.name,".chrALL.Indel",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
        ## write.table(a.indel[best.quality & an.indel,],file=paste(project.name,".chrALL.Indel_good_qual",fam[ifam],project.extension,sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
        
        
      }else{
        write.table(a.indel,file=paste(project.name,".BEST.chrALL.ACC_ACHEL",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
        #  write.table(a.indel[best.quality,],file=paste(project.name,".BEST.chrALL.ACC_ACHEL_0.025",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
        #write.table(a.indel[best.quality2,],file=paste(project.name,".BEST.chrALL.ACC_ACHEL_SUBSET",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
        # write.table(summary.reduced,file=paste(file.out.name,"ALL_GENOTYPES_analysis-maf-filtered.txt",sep="."),col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE) ## write in geno dump
        ## write.table(a.indel[an.indel,],file=paste(project.name,".chrALL.Indel",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
        ## write.table(a.indel[best.quality & an.indel,],file=paste(project.name,".chrALL.Indel_good_qual",fam[ifam],project.extension,sep=""),col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE,append=TRUE)
        
        
      }
      
    } # accumulate
    ## setwd(code.dir)
    ## source("ucsc.table.names.r")
        
    rm(a.indel)
    
    if(dont.build.summary){indels<-{}}
    ## setwd(code.dir)
    ## source("ucsc.table.names.r")
        
  } # ichr
  
} # ifam


#########################################################end of meging all the snps ##########################
##############################################################################################################
##############################################################################################################
# Separating the wanted controls to be used in the analysis.

library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(XLConnect)
library("openxlsx")
library("data.table")
library("dplyr")
library("plyr")

bad.contols <-read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/related.or.bad_IN AML.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,strip.white=TRUE,na.strings="",quote="")
all.samples.contols <-read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/sample_sheet.for.analysis_IN AML.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,strip.white=TRUE,na.strings="",quote="")
all.controls <-all.samples.contols[grepl("Control",all.samples.contols[,"SampleProject"]),"ParticipantCode"] # getting only the contol samples from the sample sheet
wanted.controls<-all.controls[!(all.controls %in% bad.contols[,"x"])] # eliminating the bad samples from the the list

#reading the merged file

code.dir <- "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/"
setwd(code.dir)
source("annotate_SNPs_subroutines.r")
core.ann <- c( "chr","start","end","ALT","REF","TYPE")
core.ann.gr <- c( "chr","start","end","strand")
#mydf  <-  read.delim("mydf_enrichment_test_file.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="")
mydf  <-  read.delim("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/2015-03-25_SKDPNexteraRun.THP.wanted.All-maf-filtered.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,strip.white=TRUE,na.strings="",quote="")
merged.snps  <-  read.csv("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/wanted_snps/merged/wanted_snps.BEST.chrALL.ACC_ACHEL.ALL.ALL_GENOTYPES_analysis.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,strip.white=TRUE,na.strings="",quote="")
#To check if the columns in merged and original mydf match or not.
#check<-cbind(key.mydf,merged.snps[,"THP-1-AS703026resist.GT"],mydf[,"THP-1-AS703026resist.GT"],merged.snps[,"THP-1-AS703026resist.AD"],mydf[,"THP-1-AS703026resist.AD"])
mydf<-merged.snps
mydf<-as.matrix(merged.snps)
a.indel<-mydf
#finaldf<-mydf  # to cbind the final coverage output
#a.indel<-as.matrix(finaldf)

#Triming leading and trailing whitespace
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
}

#a.indel <-trim(a.indel)
#a.indel[1:5,1:5]
mydf<-trim(mydf)
# chr  start     end       REF ALT
# [1,] "10" "93616"   "93616"   "C" "T"
# [2,] "10" "93694"   "93694"   "T" "C"
# [3,] "10" "1142208" "1142208" "T" "C"
# [4,] "10" "4703928" "4703928" "T" "C"
# [5,] "10" "5202104" "5202104" "C" "T"

#################

## grep("Gene.Name",a.indel[,16])
## save(list=c("column.labels"),file="column.labels.RData")
## load("column.labels.RData")
################## fast read ########### for indels **NOT for mydf!
# column.labels<-read.delim(project.file,header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
# num.vars<-dim(column.labels)[2]
# a.indel<-scan(project.file,what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
# num.lines<-length(a.indel)/(num.vars)
# dim(a.indel)<-c(num.vars,num.lines)
# a.indel<-t(a.indel)
# colnames(a.indel)<-column.labels
########################################
the.chr<-unique(a.indel[,"chr"])
print(paste("Doing Chromosome ",the.chr))

if(sum(!grepl("^chr",the.chr))>0){
  mydf[,"chr"]<-paste("chr",mydf[,"chr"],sep="")
}

#a.indel[1:5,"chr"]
#[1] "chr10" "chr10" "chr10" "chr10" "chr10"

#####Check if chr is missing 
if (grepl("^chr", mydf[,"chr"])!=TRUE){
 mydf[,"chr"]<- paste("chr",mydf[,"chr"], sep="")
}

bam.file <- "THP-1-AS703026resist.ReCal.sort.bam"
# 
# #############summary and which for mydf
all.keys <- cbind(mydf[,"chr"],mydf[,"start"],mydf[,"end"],mydf[,"strand"])
all.keys <- as.matrix(all.keys)
summary <- all.keys[c(4,6),]
summary <- all.keys
# 
# ############################################ This concatenates all the elements in the column (don't need to use in here)
# 
colnames(summary) <- core.ann.gr
#### a matrix that has position you want to look at
summary[1:5,]
# chr     start     end       strand
# [1,] "chr10" "93616"   "93616"   "+"   
# [2,] "chr10" "93694"   "93694"   "+"   
# [3,] "chr10" "1142208" "1142208" "+"   
# [4,] "chr10" "4703928" "4703928" "+"   
# [5,] "chr10" "5202104" "5202104" "+"   
# options(width=150) ## output length
# 
# data <- summary ### summary could be a subset of a.indel - . a.indel[wanted,core.ann.gr)
# data.gr <- GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])
# 
# which <-   data.gr
# which
# GRanges object with 10438 ranges and 0 metadata columns:
#   seqnames               ranges strand
# <Rle>            <IRanges>  <Rle>
#   [1]    chr10   [  93616,   93616]      +
#   [2]    chr10   [  93694,   93694]      +
#   [3]    chr10   [1142208, 1142208]      +
#   [4]    chr10   [4703928, 4703928]      +
#   [5]    chr10   [5202104, 5202104]      +
#   ...      ...                  ...    ...
# [10434]     chrY [15448218, 15448218]      +
#   [10435]     chrY [15448218, 15448218]      +
#   [10436]     chrY [21894635, 21894636]      +
#   [10437]     chrY [21894635, 21894636]      +
#   [10438]     chrY [21894636, 21894636]      +
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths


#############summary and which for merged.snps
# all.keys <- cbind(merged.snps[,"chr"],merged.snps[,"start"],merged.snps[,"end"],merged.snps[,"strand"])
# all.keys <- as.matrix(all.keys)
# summary <- all.keys[c(4,6),]
# summary <- all.keys
# summary
# summary[1:5,]
# chr     start     end       strand
# [1,] "chr10" "93616"   "93616"   "+"   
# [2,] "chr10" "93694"   "93694"   "+"   
# [3,] "chr10" "1142208" "1142208" "+"   
# [4,] "chr10" "4703928" "4703928" "+"   
# [5,] "chr10" "5202104" "5202104" "+"   
############################################ This concatenates all the elements in the column (don't need to use in here)

# colnames(summary) <- core.ann.gr
# summary #### a matrix that has position you want to look at
# options(width=150) ## output length

data <- summary ### summary could be a subset of a.indel - . a.indel[wanted,core.ann.gr)
data.gr <- GRanges(seqnames =data[,"chr"],ranges = IRanges(start=as.numeric(data[,"start"]),end=as.numeric(data[,"end"])),strand=data[,"strand"])

which <-   data.gr
which
######
#bam.dir <- "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AMAS/BAM"
#bam.dir <- c("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AMAS/BAM","/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/MODY/BAM/","/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AOGC-NGS/BAM/")
#bam.dir <-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/MODY/BAM/"
#bam.dir
#my.dir<-setwd("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/")
#setwd("./AMAS/BAM")
#all.bam.files<-list.files(getwd())
#setwd('../')
#mydf <- within(mydf, start  <-  paste(chr,start, sep=':'))
mydf[,"start"]<-paste(mydf[,"chr"],mydf[,"start"],sep=":")
#mydf<-as.matrix(mydf)     ##########format of the a.indel or this is how a.indel should be
mydf[1:10,1:9]  #This is a.indel
# chr     start           end         REF ALT TYPE  refGene::location   refGene::type                                 refGene::gene
# 1  "chr10" "chr10:93616"   "    93616" "C" "T" "snp" "nonsynonymous SNV" "TUBB8:NM_177987:exon4:c.G716A:p.C239Y,"      "TUBB8"      
# 2  "chr10" "chr10:93694"   "    93694" "T" "C" "snp" "nonsynonymous SNV" "TUBB8:NM_177987:exon4:c.A638G:p.K213R,"      "TUBB8"      
# 3  "chr10" "chr10:1142208" "  1142208" "T" "C" "snp" "intronic"          "WDR37"                                       "WDR37"      
# 4  "chr10" "chr10:4703928" "  4703928" "T" "C" "snp" "ncRNA_exonic"      "LINC00705"                                   "LINC00705"  
# 5  "chr10" "chr10:5202104" "  5202104" "C" "T" "snp" "ncRNA_exonic"      "AKR1CL1"                                     "AKR1CL1"    
# 6  "chr10" "chr10:5255025" "  5255025" "A" "G" "snp" "nonsynonymous SNV" "AKR1C4:NM_001818:exon7:c.A749G:p.Q250R,"     "AKR1C4"     
# 7  "chr10" "chr10:5683885" "  5683885" "A" "C" "snp" "nonsynonymous SNV" "ASB13:NM_024701:exon5:c.T557G:p.L186R,"      "ASB13"      
# 8  "chr10" "chr10:5790420" "  5790420" "T" "C" "snp" "nonsynonymous SNV" "FAM208B:NM_017782:exon15:c.T5036C:p.V1679A," "FAM208B"    
# 9  "chr10" "chr10:7601836" "  7601836" "C" "T" "snp" "UTR3"              "ITIH5"                                       "ITIH5"      
# 10 "chr10" "chr10:7763661" "  7763661" "A" "G" "snp" "nonsynonymous SNV" "ITIH2:NM_002216:exon8:c.A788G:p.N263S,"      "ITIH2" 

#for (j in 1:length(bam.dir)){
#all.bam.files<-list.files(bam.dir[j])
#all.bam.names<-sub(".*//","",all.bam.files)
#all.bam.names<-all.bam.names[grepl(".bam$",all.bam.names)]
#all.bam.files<-all.bam.files[grepl(".bam$",all.bam.files)]
#all.bam.files[1:10]
# [1] "CHI_004.2_C4UAMACXX.ReCal.sort.bam"      "CHI_007.3_C6AR7ACXX-1-16.ReCal.sort.bam" "CHI_012.3_C4UAMACXX.ReCal.sort.bam"     
# [4] "CHI_017.1.ReCal.sort.bam"                "CHI_017.2.ReCal.sort.bam"                "CHI_024.3_C4UAMACXX.ReCal.sort.bam"     
# [7] "CHI_035.3_C4UAMACXX.ReCal.sort.bam"      "CHI_042.1_C6AR7ACXX-1-14.ReCal.sort.bam" "CHI_042.2_C6AR7ACXX-1-15.ReCal.sort.bam"
# [10] "CHI_042.3_C4UAMACXX.ReCal.sort.bam"     

#all.bam.names<-all.bam.names[grepl(".bam$",all.bam.names)&grepl("^THP",all.bam.files)]
#all.wanted.bam<-c(wanted.controls,"THP-1-AS703026resist.ReCal.sort.bam","THP-1-normal-parent.ReCal.sort.bam","THP-1-Selumetinibresist.ReCal.sort.bam") #selecting the contols and the THP BAM files 
#all.bam.files <- list.files(bam.dir,full.names = TRUE, recursive = TRUE)
#wanted.names <- paste0(all.wanted.bam, collapse="|")
#all.bam.files.with.control<-grep(wanted.names, all.bam.files, value=TRUE)
#all.bam.files.with.control<-all.bam.files.with.control[grepl(".bam$",all.bam.files.with.control)]
#unlist(mapply(function(x, y) grep(y, x, value = TRUE), x = all.bam.files, y = wanted.controls))

#controlA<-paste("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AMAS/BAM//",wanted.controls,sep="")
#controlB<-paste("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/MODY/BAM//",wanted.controls,sep="")
#controlC<-paste("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AOGC-NGS/BAM/",wanted.controls,sep="")
#control.paths.all<-c(controlA,controlB,controlC)
#control.paths.all<-paste(control.paths.all,wanted.controls,sep="")
#grepl(/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/MODY/BAM//,all.bam.files
#all.bam.files<-all.bam.files[grepl(".bam",all.bam.files)]
#all.bam.files<-all.bam.files[grepl(".bam",all.bam.files)&grepl("^THP",all.bam.files)]
#all.bam.filesa<-all.bam.files[match(wanted.controls[136:236],all.bam.files)&grepl("^THP",all.bam.files)]
#list.files(list.files(".",pattern="/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AMAS/BAM"|"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/MODY/BAM/"|"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AOGC-NGS/BAM/"))
##########################
#all.files<-Vectorize(list.files)(bam.dir)
#unlisted<-unlist(all.files)
#unlisted.bam<-unlisted[grepl(".bam$",unlisted)]
#51623174

#all.bam.files<-all.bam.files.with.control  #all.bam.files contains all the bam files to be looped over for the samples and controls

params <- ScanBamParam(which=which,flag=scanBamFlag(isUnmappedQuery=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE),simpleCigar = FALSE,reverseComplement = FALSE,what=c("qname","flag","rname","seq","strand","pos","qwidth","cigar","qual","mapq") )  ### NOTE isValidVendorRead=FALSE shoudl be TRUE


param.pile <- PileupParam(max_depth=2500, min_base_quality=20, min_mapq=0,min_nucleotide_depth=1, min_minor_allele_depth=0,distinguish_strands=TRUE, distinguish_nucleotides=TRUE,ignore_query_Ns=TRUE, include_deletions=TRUE,cycle_bins=numeric() )
#save.merged<-merged.snps
#merged.snps<-(within(merged.snps, start  <-  paste(chr,start, sep=':')))
#merged.snps<-as.matrix(merged.snps)
#test<-merged.snps
#trim <- function( x ) {
#  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x)
#}
#test<-trim(test)
#bam.file<-all.bam.files[i]
#bam.file<-c("THP-1-Selumetinibresist.ReCal.sort.bam")

#setwd(bam.dir)
#############################################################################
## Function that does all the work
#bam.file<-all.bam.files[1] #to check if the function works or not
#setwd(bam.dir)
func.bamAD <- function (bam.file,mydf){

test <- pileup(bam.file,scanBamParam=params,pileupParam=param.pile,ApplyPileupsParam=param)
#setwd("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/")

newtest<-cbind(start=paste(test[,"seqnames"],test[,"pos"],sep=":"),test)
#newtest <- within(test, start  <-  paste(seqnames,pos, sep=':'))
newtest[1:5,]

yyy <- tapply(newtest[,"count"],list(newtest[,"nucleotide"],newtest[,"which_label"]),sum)
dim(yyy)
yyy[1:5,]

coverage.bam <- t(yyy)
coverage.bam <- cbind(rownames(coverage.bam),coverage.bam)
colnames(coverage.bam)[1] <- c("start")
rownames(coverage.bam) <- NULL
coverage.bam[1:10]
# [1] "chr10:93616-93616.1"      "chr10:93694-93694.2"      "chr10:1142208-1142208.3"  "chr10:4703928-4703928.4"  "chr10:5202104-5202104.5" 
# [6] "chr10:5255025-5255025.6"  "chr10:5683885-5683885.7"  "chr10:5790420-5790420.8"  "chr10:7601836-7601836.9"  "chr10:7763661-7763661.10"

#coverage.bam<-(coverage.bam,c(gsub("*-.*", "\\1", coverage.bam[,"which_label"]),"newcol"))
coverage.bam[,"start"] <-c(gsub("*-.*", "\\1", coverage.bam[,"start"])) #removes the numbers after "-"

#Since the bam file is missing some of the ncRNA_intronic, UTR, need to remove the difference between the file from the coverage.bam matrix  


###########

dim(mydf)  #6179,268
dim(coverage.bam) # 6017,8
missing <- !mydf[,"start"]%in%coverage.bam[,"start"]
sum(missing) # 148(excluding the duplicates)

rows.present <- mydf[(mydf[,"start"]%in%coverage.bam[,"start"]),] #matching rows including duplicates (i.e only -148, still 14X2 are duplicates)
dim(rows.present)

dim(coverage.bam)


coverage.bam <- coverage.bam[,c("start","A","C","G","T","N","=","-") ]
coverage.bam[1:5,]
# start           A  C    G  T    N  =  - 
#   [1,] "chr10:93616"   NA "65" NA "13" NA NA NA
# [2,] "chr10:93694"   NA "24" NA "57" NA NA NA
# [3,] "chr10:1142208" NA "17" NA NA   NA NA NA
# [4,] "chr10:4703928" NA NA   NA NA   NA NA NA
# [5,] "chr10:5202104" NA NA   NA "17" NA NA NA
#######################
result<-cbind(coverage.bam,rows.present)
result<-result[!duplicated(result),]
cols <- colnames(result) == "start"
colnames(result)[cols] <- paste0("start", seq.int(sum(cols)))
#check if they are same
colnames(result)[which(colnames(result) == "start1")]<-"start"

sum(result[,"start"] %in% result[,"start2"])
wanted.pos<-!grepl("^indel",result[,"TYPE"]) # removing all the indel lines
sum(wanted.pos)
wanted.result<-result[wanted.pos,]
#newdf<-bamAD(wanted.result)
#wanted.result[,"start2"]<- NULL
want.result<-wanted.result  ##has to be matrix to use the function newBamAD
df <- na.omit(want.result)
types <- c("A", "T", "G", "C")
df[types][is.na(df[types])] <- 0
########function to get ADs 
newBamAD <- function (x,base.types=c("A","C","G","T")) {
  # the version above
  rownames(x) <- 1:nrow(x)
  ref <- x[cbind(1:nrow(x), x[, 'REF'])]
  alt <- x[cbind(1:nrow(x), x[, 'ALT'])]
  which.flat <- grep('flat$', x[, 'TYPE'])
  
  alt[which.flat] <-  sapply(which.flat, function (i,base.types) {
   sum(as.numeric(x[i, c( base.types[!( base.types %in% x[i, 'REF'])] )] ) ,na.rm=TRUE) },base.types)
  cbind(x[,c("start","REF","ALT","TYPE")],bam.AD=paste(ref, alt, sep=','))
 # cbind(x, bam.AD=paste(ref, alt, sep=','))
}

##############################

myvec<-newBamAD(want.result)
#colnames(myvec)[5]<-paste(bam.file,"FAD",sep=".")  ##change the column bam.AD to all file name of respective bam file.
colnames(myvec)[5]<-gsub(".*BAM/", "", bam.file)
#colnames(myvec)[5]<-paste(unlist(strsplit(bam.file, "_[^_]+$|.ReCal"))[1],"FAD",sep=".")
return(myvec)
}


#paste(unlist(strsplit(all.bam.files[22], "_[^_]+$"))[1],"FAD",sep=".")
#paste(unlist(strsplit(all.bam.files[22], "_C|.ReCal"))[1],"FAD",sep=".")
#unlist(strsplit(all.bam.files[1], "_[A-Za-z0-9]+$"))

bam.dir <- c("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AMAS/BAM",
             "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/MODY/BAM/",
             "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/AOGC-NGS/BAM/",
             "/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Projects/SDDS/BAM/")
bam.dir

#all.wanted.bam<-c(wanted.controls,"THP-1-AS703026resist","THP-1-normal-parent","THP-1-Selumetinibresist")



#wanted.names <- paste0(all.wanted.bam, collapse="|")
#missing<-is.na(pmatch(wanted.controls,all.bam.files,nomatch = NA_integer_, duplicates.ok = FALSE))
#wanted.controls[missing]


# all.bam.files<-{}
# not.found<-{}
# found<-{}
# 
# #j<-2
# for (j in 1:length(bam.dir)){
#   #print(paste("Processing Bam directory:",bam.dir[j]))
#   print(paste(paste("Directory number",(j),"Processing Bam directory:",sep=":"),bam.dir[j]))
#   all.bam.files<-list.files(bam.dir[j])
#   all.bam.files<-all.bam.files[grepl(".bam$",all.bam.files)]
#   all.bam.files<-c(all.bam.files,unlist(lapply(all.bam.files, function(bam) all.wanted.bam[grepl(bam, all.wanted.bam)])))
#   #track.record <- unlist(as.matrix(lapply(all.wanted.bam, function(bam) grep(bam, all.bam.files)))==0)
#   #indx2 <- unlist(is.na(track.record[,1]))
#   indx2<- unlist(lapply(lapply(all.wanted.bam, function(bam) agrepl(bam, all.bam.files)),any))
#   not.found<-all.wanted.bam[indx2==TRUE]
#   found<-all.wanted.bam[indx2==FALSE]
#   #wanted.file.names<-all.wanted.bam[pmatch(all.wanted.bam,all.bam.files)]
#   #missing<-is.na(wanted.file.names)
#   #wanted.file.names<-wanted.file.names[!missing]
#   #all.file.names<-c(all.file.names,wanted.file.names)
#   #all.file.names<-c(all.file.names,grep(wanted.names, all.bam.files, value=TRUE))
#   #indx <- unlist(lapply(all.wanted.bam, function(bam) grep(bam, all.bam.files)))
#   #print (indx)
#   #found<-c(found,all.wanted.bam[indx])
#   #not.found<-c(not.found,all.wanted.bam[-indx])
# }



#Since all the samples do not have .GT columns, need to weed out those without .GT columns from all.samples
wanted.columns.indel.GT<-colnames(a.indel)[grepl(".GT$",colnames(a.indel))]
all.wanted.samples<-c(wanted.controls,"THP-1-AS703026resist","THP-1-normal-parent","THP-1-Selumetinibresist") 
###selecting the contols and the THP BAM files, wanted.controls have been worked out above 
all.wanted.samples.GT<-paste(all.wanted.samples,"GT",sep=".")
#all.wanted.samples.GT[(is.na(match(all.wanted.samples.GT ,wanted.columns.indel)))]
##finding the common samples in wanted.columns.indel.GT and all.wanted.samples.GT
#intersect(wanted.columns.indel.GT,all.wanted.samples.GT)
wanted.samples.GT<-wanted.columns.indel.GT[(wanted.columns.indel.GT%in%all.wanted.samples.GT)]
wanted.samples.FAD<-gsub(".GT",".FAD",wanted.samples.GT)
wanted.samples<-gsub(".GT","",wanted.samples.GT)
length(wanted.samples)

#total.FAD is where all.FAD are column bind and becomes a large matrix
total.FAD<-{}
for (j in 1:length(bam.dir)){
  #print(paste("Processing Bam directory:",bam.dir[j]))
  print(paste(paste("Directory number",(j),"Processing Bam directory:",sep=":"),bam.dir[j]))
  all.bam.files<-list.files(bam.dir[j])
  #all.bam.names<-sub(".*//","",all.bam.files)
  #all.bam.names<-all.bam.names[grepl(".bam$",all.bam.names)]
  all.bam.files<-all.bam.files[grepl(".bam$",all.bam.files)]
  print(paste("*************The total number of bam files in this dir:", (length(all.bam.files)),"*************"))
  all.bam.files[1:10]
  # [1] "CHI_004.2_C4UAMACXX.ReCal.sort.bam"      "CHI_007.3_C6AR7ACXX-1-16.ReCal.sort.bam" "CHI_012.3_C4UAMACXX.ReCal.sort.bam"     
  # [4] "CHI_017.1.ReCal.sort.bam"                "CHI_017.2.ReCal.sort.bam"                "CHI_024.3_C4UAMACXX.ReCal.sort.bam"     
  # [7] "CHI_035.3_C4UAMACXX.ReCal.sort.bam"      "CHI_042.1_C6AR7ACXX-1-14.ReCal.sort.bam" "CHI_042.2_C6AR7ACXX-1-15.ReCal.sort.bam"
  # [10] "CHI_042.3_C4UAMACXX.ReCal.sort.bam"     
 
  wanted.names <- paste0(wanted.samples, collapse="|")
  all.bam.files.with.control<-grep(wanted.names, all.bam.files, value=TRUE)
  all.bam.files<-all.bam.files.with.control
 print(paste("The wanted number of bam files in this directory:", (length(all.bam.files))))
 if(length(all.bam.files)==0){
   next
 }else{
 setwd(bam.dir[j])
}
all.FAD<- {}
  for(i in 1:length(all.bam.files)){
  #print(paste((1)"Doing bam file:",all.bam.files[i]))
  print(paste(paste("File number",(i),"Doing bam file:",sep=":"),all.bam.files[i]))
  an.AD = func.bamAD(all.bam.files[i],mydf)
  if(is.null(all.FAD)){all.FAD<-an.AD
    }else{
  all.FAD<-cbind(all.FAD,an.AD)
    }
  print(paste("The dimension of the all.FAD table is:",dim(all.FAD)[1],"rows","and",dim(all.FAD)[2],"columns"))
  }
if(is.null(total.FAD)){total.FAD<-all.FAD
  }else{
total.FAD<-cbind(total.FAD,all.FAD)
  }
print(paste("###########The dimension of the total.FAD table is:",dim(total.FAD)[1],"and",dim(total.FAD)[2],"columns","###########"))
}

setwd("/media/TRI-T-DRIVE-taneupan/uq/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/")
write.table(total.FAD,file="/media/TRI-T-DRIVE-taneupan/uq/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/totalFAD_Aug16.txt",col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE, na="")


GG<-total.FAD  #saving the all.FAD in case needed again!
total.FAD<-GG

#bam.colms<-colnames(total.FAD)[grepl(".FAD$",colnames(total.FAD))]
#This is the wanted FAD file
bam.colms<-grep(wanted.names, colnames(total.FAD), value=TRUE)


bamcolumns <- unique(sub("(.*)(_.*)", "\\1", bam.colms))
bamcolumns <- bamcolumns[bamcolumns%in%wanted.samples][c(1:48,54:73)]
wantedcols <- paste0(bamcolumns, collapse="|")
bam.colms<-grep(wantedcols, colnames(total.FAD), value=TRUE)
thp.samples<-c("THP-1-normal-parent.ReCal.sort.bam","THP-1-AS703026resist.ReCal.sort.bam","THP-1-Selumetinibresist.ReCal.sort.bam")
bam.colms<-c(thp.samples,bam.colms)
bam.FAD<-total.FAD[,c("start","REF","ALT","TYPE",bam.colms)]

# bam.FAD<-gsub(',NA',',0',bam.FAD)
# bam.FAD<-gsub('NA,','0,',bam.FAD)


snp.only<-grepl("^snp",a.indel[,"TYPE"])
wanted.indel<-a.indel[snp.only,]

key.wanted.indel <- paste(paste("chr",wanted.indel[,"chr"],sep=""), wanted.indel[,"start"],wanted.indel[,"REF"], wanted.indel[,"ALT"], wanted.indel[,"TYPE"], sep=":")
key.wanted.indel[1:10]
# [1] "chr10:93616:C:T:snp"   "chr10:93694:T:C:snp"   "chr10:1142208:T:C:snp" "chr10:4703928:T:C:snp" "chr10:5202104:C:T:snp" "chr10:5255025:A:G:snp" "chr10:5683885:A:C:snp"
# [8] "chr10:5790420:T:C:snp" "chr10:7601836:C:T:snp" "chr10:7763661:A:G:snp"
#check if they align correctly to be cbind
key.FAD <- paste(bam.FAD[,"start"], bam.FAD[,"REF"],all.FAD[,"ALT"],bam.FAD[,"TYPE"],sep=":")
key.FAD[1:10]
# [1] "chr10:93616:C:T:snp"   "chr10:93694:T:C:snp"   "chr10:1142208:T:C:snp" "chr10:4703928:T:C:snp" "chr10:5202104:C:T:snp" "chr10:5255025:A:G:snp" "chr10:5683885:A:C:snp"
# [8] "chr10:5790420:T:C:snp" "chr10:7601836:C:T:snp" "chr10:7763661:A:G:snp"
posns.wanted.indel <- match(key.wanted.indel,key.FAD)
sum(is.na(posns.wanted.indel)) #this has to be zero
############################################################################################################################################j###########
bam.FAD<-cbind('refGene::type'=wanted.indel[,"refGene::type"],bam.FAD)
#key.FAD[1:10]
# [1] "chr10:93616:C:T:snp:TUBB8:NM_177987:exon4:c.G716A:p.C239Y,"        "chr10:93694:T:C:snp:TUBB8:NM_177987:exon4:c.A638G:p.K213R,"       
# [3] "chr10:1142208:T:C:snp:WDR37"                                       "chr10:4703928:T:C:snp:LINC00705"                                  
# [5] "chr10:5202104:C:T:snp:AKR1CL1"                                     "chr10:5255025:A:G:snp:AKR1C4:NM_001818:exon7:c.A749G:p.Q250R,"    
# [7] "chr10:5683885:A:C:snp:ASB13:NM_024701:exon5:c.T557G:p.L186R,"      "chr10:5790420:T:C:snp:FAM208B:NM_017782:exon15:c.T5036C:p.V1679A,"
# [9] "chr10:7601836:C:T:snp:ITIH5"                                       "chr10:7763661:A:G:snp:ITIH2:NM_002216:exon8:c.A788G:p.N263S,"   

#key.indel<-paste(a.indel[,"chr"],a.indel[,"start"],a.indel[,"REF"],a.indel[,"ALT"],a.indel[,"TYPE"],a.indel[,"refGene::type"], sep=":")
#key.indel[1:10]
# [1] "chr10:93616:C:T:snp:TUBB8:NM_177987:exon4:c.G716A:p.C239Y,"        "chr10:93694:T:C:snp:TUBB8:NM_177987:exon4:c.A638G:p.K213R,"       
# [3] "chr10:1142208:T:C:snp:WDR37"                                       "chr10:4703928:T:C:snp:LINC00705"                                  
# [5] "chr10:5202104:C:T:snp:AKR1CL1"                                     "chr10:5255025:A:G:snp:AKR1C4:NM_001818:exon7:c.A749G:p.Q250R,"    
# [7] "chr10:5683885:A:C:snp:ASB13:NM_024701:exon5:c.T557G:p.L186R,"      "chr10:5790420:T:C:snp:FAM208B:NM_017782:exon15:c.T5036C:p.V1679A,"
# [9] "chr10:7601836:C:T:snp:ITIH5"                                       "chr10:7763661:A:G:snp:ITIH2:NM_002216:exon8:c.A788G:p.N263S," 

#posns<-match(key.indel,key.FAD) if to be mapped with the original indel
key.FAD <-paste(bam.FAD[,"start"], bam.FAD[,"REF"],all.FAD[,"ALT"],bam.FAD[,"TYPE"],bam.FAD[,"refGene::type"],sep=":")
key.FAD[1:10]
# [1] "chr10:93616:C:T:snp:TUBB8:NM_177987:exon4:c.G716A:p.C239Y,"        "chr10:93694:T:C:snp:TUBB8:NM_177987:exon4:c.A638G:p.K213R,"       
# [3] "chr10:1142208:T:C:snp:WDR37"                                       "chr10:4703928:T:C:snp:LINC00705"                                  
# [5] "chr10:5202104:C:T:snp:AKR1CL1"                                     "chr10:5255025:A:G:snp:AKR1C4:NM_001818:exon7:c.A749G:p.Q250R,"    
# [7] "chr10:5683885:A:C:snp:ASB13:NM_024701:exon5:c.T557G:p.L186R,"      "chr10:5790420:T:C:snp:FAM208B:NM_017782:exon15:c.T5036C:p.V1679A,"
# [9] "chr10:7601836:C:T:snp:ITIH5"                                       "chr10:7763661:A:G:snp:ITIH2:NM_002216:exon8:c.A788G:p.N263S,"     
key.wanted.indel <- paste(paste("chr",wanted.indel[,"chr"],sep=""),wanted.indel[,"start"],wanted.indel[,"REF"],wanted.indel[,"ALT"],wanted.indel[,"TYPE"],wanted.indel[,"refGene::type"],sep=":")
key.wanted.indel[1:10]
# [1] "chr10:93616:C:T:snp:TUBB8:NM_177987:exon4:c.G716A:p.C239Y,"        "chr10:93694:T:C:snp:TUBB8:NM_177987:exon4:c.A638G:p.K213R,"       
# [3] "chr10:1142208:T:C:snp:WDR37"                                       "chr10:4703928:T:C:snp:LINC00705"                                  
# [5] "chr10:5202104:C:T:snp:AKR1CL1"                                     "chr10:5255025:A:G:snp:AKR1C4:NM_001818:exon7:c.A749G:p.Q250R,"    
# [7] "chr10:5683885:A:C:snp:ASB13:NM_024701:exon5:c.T557G:p.L186R,"      "chr10:5790420:T:C:snp:FAM208B:NM_017782:exon15:c.T5036C:p.V1679A,"
# [9] "chr10:7601836:C:T:snp:ITIH5"                                       "chr10:7763661:A:G:snp:ITIH2:NM_002216:exon8:c.A788G:p.N263S,"     
posns<-match(key.wanted.indel,key.FAD)
missing<-is.na(posns)
sum(missing) #has to be zero
#test<-bam.FAD
####This loop replaces the column names in the all.FAD matrix with the names in all.wanted.bam
mycolumns<-colnames(bam.FAD)
##this loop generated foul result while replacing column names
#for(i in as.character(all.wanted.bam)){
#  mycolumns[grepl(i, mycolumns)] <- i
#}
#Replacing column names
library(gsubfn)
f <- function(x,y,z) if (z=="_") y else strsplit(x, "_", fixed=T)[[1]][[1]]
new.bam.cols<-gsubfn("([^_]+_[^_]+)(.).*", f, mycolumns, backref=2)
new.bam.cols<-gsub(".ReCal.sort.bam" ,"",new.bam.cols)
#These three lines are worked out to get controls and samples for next analysis of pvalue enrichment
wanted.samples<-new.bam.cols[c(6:66)] #this is going to be used for pvalue below
new.bam.cols<-c(new.bam.cols[1:5],paste(new.bam.cols[-c(1:5)],"FAD",sep=".")) #back to the replacement of column names
colnames(bam.FAD)<-new.bam.cols
all.samples.in.all.FAD<-colnames(bam.FAD)[!colnames(bam.FAD)%in%c("refGene::type","start","REF","ALT","TYPE")]
copy.bam.FAD<-cbind(key.FAD=paste(bam.FAD[,"start"],bam.FAD[,"REF"],bam.FAD[,"ALT"],bam.FAD[,"TYPE"],bam.FAD[,"refGene::type"],sep=":"),bam.FAD)
#check<-colnames(all.FAD)[!colnames(all.FAD)%in%c("refGene::type","start","REF","ALT","TYPE")]
#verify.if.correct<-paste(all.samples.in.all.FAD,check,sep="@@@@@")
#tru<-all.wanted.bam%in%all.samples.in.all.FAD
#verify<-paste(tru,verify.if.correct,sep=":") #check if the samples are missing as I have not included all the directories to loop over all.wanted.bam
#after verifying the missing samples
#included.column<-all.samples.in.all.FAD[c(1:16,22:91,106:127,197:218)] #Manually Selecting the correct names (ones without ambiguity)

#write.table(verify,file="/media/TRI-T-DRIVE-taneupan/uq/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/verify.txt",col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE, na="")

#wanted.FAD<-copy.test[,c("startFAD",included.column)]
#colnames(wanted.FAD)<-paste(colnames(wanted.FAD),"FAD",sep=".")
#bam.columns.AD <-paste(included.column,"FAD",sep=".")
#a.indel<-cbind(a.indel,bam.FAD[posns,bam.columns.AD])
wanted.indel<-cbind(key.indel=paste(wanted.indel[,"chr"],wanted.indel[,"start"],wanted.indel[,"REF"],wanted.indel[,"ALT"],wanted.indel[,"TYPE"],wanted.indel[,"refGene::type"],sep=":"),wanted.indel)
sum(!key.wanted.indel%in%key.FAD) # the mismatch between key.wanted.indel and key.FAD should be zero 
#To check the AD columns with original AD columns (subsetting the matrix)
#check<-final.indel[,c(1:14,grep(".AD$",colnames(final.indel)))]
final.FAD<-cbind(key.start=copy.bam.FAD[,"key.FAD"],wanted.indel,copy.bam.FAD[,all.samples.in.all.FAD])
#final.FAD[c(1:10,1123:1128]
write.table(final.FAD,file="/media/TRI-T-DRIVE-taneupan/uq/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/final.FAD.txt",col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE, na="")


#To make sure we have a .GT column for every .FAD cols
check.colmns<-bam.colms
gtcolmns<-colnames(final.FAD)[grepl(".GT$",colnames(final.FAD))]
check.colmns

f <- function(x,y,z) if (z=="_") y else strsplit(x, "_", fixed=T)[[1]][[1]]
check.colmns<-gsubfn("([^_]+_[^_]+)(.).*", f, check.colmns, backref=2)
check.colmns<-gsub(".ReCal.sort.bam" ,"",check.colmns)
check.colmns<-paste(check.colmns,"GT",sep=".")
(check.colmns%in%gtcolmns) # has to be all true

###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################
#######################################################   END of Pile up !!! Begins get_pvalue_for_enrichment  #################################################################
###############################################################################################################################################
###############################################################################################################################################
###############################################################################################################################################

#load("/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-03-16_AllAMLandLung/Analysis/FAD.recovered.RData")

#all.samples.contols
#### NEED GT AND FAD IN THE SAME MATRIX
## genotypes<-a.indel[,paste(all.possible.samples,"GT",sep=".")] #### use a.indel.ori otherwise
## a.indel.stats<-cbind(a.indel.stats,genotypes)
#pheno.ori[1:10,1:10]
# SampleProject FamilyCode SAMPLE PaternalID MaternalID Sex AffectionStatus UPN Sample.Number sample.Source
# 32             1        ALL      1          0          0  -9               2   1             1           Tom
# 33             1        ALL     10          0          0  -9               2   9             9           Tom
# 34             1        ALL    100          0          0  -9               2  96            96           Tom
# 35             1        ALL    101          0          0  -9               2  97            97           Tom
# 37             0        ALL  10694          0          0  -9               1  NA            NA          <NA>
#   38             0        ALL  10695          0          0  -9               1  NA            NA          <NA>
#   40             0        ALL  11180          0          0  -9               1  NA            NA          <NA>
#   41             0        ALL  11181          0          0  -9               1  NA            NA          <NA>
#   43             0        ALL  11286          0          0  -9               1  NA            NA          <NA>
#   44             0        ALL  11287          0          0  -9               1  NA            NA          <NA>
#**********************************************************************
#cellularity<- pheno.ori[,c("SAMPLE","Blast.Count..")]
##cellularity<-sample.types[,c("ParticipantCode","Blast.Count..")] #AN-my samples
# > pheno.ori[,c("SAMPLE","Blast.Count..")][1:10,]
# SAMPLE Blast.Count..
# 32      1            76
# 33     10            55
# 34    100            92
# 35    101            91
# 37  10694            NA
# 38  10695            NA
# 40  11180            NA
# 41  11181            NA
# 43  11286            NA
# 44  11287            NA
#cellularity[,"Blast.Count.."]<-cellularity[,"Blast.Count.."]/100
#cellularity[is.na(cellularity[,"Blast.Count.."]) & pheno.ori[,"Control"] ,"Blast.Count.."]<-1 ## unnoen then make pure
#**********************************************************************
#cellularity<-all.samples.contols[,c("SAMPLE","Blast.Count..")]
####
#These are all the samples (both controls and cancer samples to be included in Pvalue analysis below)
cancer<-wanted.samples[grepl("^THP",wanted.samples)] #These are the cancer samples to be used in Pvalue analysis below
#Now checking if the .GT and >FAD columns are present for all those cancer and control samples

cancer.GT<-paste(cancer,"GT",sep=".")
Controls<-wanted.samples[!grepl("^THP",wanted.samples)] #These are the Controls to be used in Pvalue analysis below
length(Controls)
Controls.GT<-paste(Controls,"GT",sep=".")

wanted.samplesFAD<-paste(wanted.samples,"FAD",sep=".")
wanted.samplesGT<-paste(wanted.samples,"GT",sep=".")
wanted.samplesAD<-paste(wanted.samples,"AD",sep=".")
#indels<-final.FAD[,c(colnames(final.FAD)[1:17],wanted.samplesGT,wanted.samplesFAD)]
indels<-final.FAD

indels[] <- lapply(indels, sub, pattern = "NA(?=,[0-9])|(?<=[0-9],)NA", 
                      replacement = "0", perl = TRUE)


#cellularity<-read.delim("cellularity.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)
cellularity<-cbind(Samples=wanted.samples, Cellularity=1)
cellularity[1:5,]


# Sample Cellularity
# 1       THP-1-AS703026resist           1
# 2    THP-1-Selumetinibresist           1

## 10  02851          1.00
## 14  03199          1.00
## 21  03259          1.00
## 25  03345          1.00
## 32      1          0.76

#############################
#cellularity<- pheno.ori[,c("SAMPLE","Blast.Count..")]
#cellularity[,"Blast.Count.."]<-cellularity[,"Blast.Count.."]/100
#cellularity[is.na(cellularity[,"Blast.Count.."]) & pheno.ori[,"Control"] ,"Blast.Count.."]<-1 ## unnoen then make pure
# cellularity[pheno.ori[,"AML"],]


PDs<-pheno.ori[pheno.ori[,"AML-Child"] | pheno.ori[,"Asian-AML-Child"] | pheno.ori[,"Asian-AML"]  | pheno.ori[,"AML-NotDiagnosis-Child"] | pheno.ori[, "Asian-AML-NotDiagnosis-Child"],"SAMPLE"]
## PDs<-c(PDs,"LPH-001-27_PD")   table(pheno.ori[pheno.ori[,"AffectionStatus"]==2,"SAMPLE"] )
PDs ## other samplesyou want to correct



#cancer<-pheno.ori[pheno.ori[,"AML"],"SAMPLE"]
#cancer
## [131] "AMLM12039M-A"    "AMLM12040A-J"    "AMLM12PAH016K-B"
#******
cancer
# "THP-1-normal-parent"     "THP-1-AS703026resist"    "THP-1-Selumetinibresist"

#cancer.alt.counts<-alt.reads.reference.calls(a.indel,cancer,threshold=1)
#cancer.alt.counts.true<-alt.reads.Non.reference.calls(a.indel,cancer,threshold=1)

#Controls<-pheno.ori[pheno.ori[,"Control"],"SAMPLE"]
## > Controls
##   [1] "02851"          "03199"          "03259"          "03345"          "10694"          "10695"          "11180"          "11181"          "11286"          "11287"          "13371"         
##  [12] "13372"          "13374"          "13375"          "13377"          "13378"          "13395"          "13396"          "13398"          "13399"          "13401"          "13402"         
##  [23] "13425"          "13426"          "13428"          "13431"          "13432"          "13436"          "13437"          "13439"          "13440"          "13479"          "13480"         
##  [34] "13482"          "13483"          "13486"          "13487"          "3170"           "3212"           "3279"           "3291"           "3294"           "3295"           "3309"          
##  [45] "3310"           "34.8"   

#********
Controls[1:10]
# [1] "10694" "10695" "11180" "11181" "11286" "11287" "13371" "13372" "13374" "13375"         

a.indel.stats[1:5,1:20]
#Control.alt.counts<-alt.reads.reference.calls(a.indel.ori,Controls,AD.extension="AD",threshold=1)
## > a.indel.stats[1:5,1:20]
##                                                           id chr    start      end REF ALT TYPE  GENE 02851.TAD 02851.FAD 02851.DUP 03199.TAD 03199.FAD 03199.DUP 03259.TAD 03259.FAD 03259.DUP
## chr10:50942493:50942493:C:T:snp 10:50942493:50942493:C:T:snp  10 50942493 50942493   C   T  snp OGDHL       2,0       1,0       0,0       0,0       0,0       0,0       2,0       2,0       0,0
## chr10:50942558:50942558:G:A:snp 10:50942558:50942558:G:A:snp  10 50942558 50942558   G   A  snp OGDHL       1,0       1,0       0,0       0,2       0,1       1,1       1,0       1,0       0,0
## chr10:50942763:50942763:C:T:snp 10:50942763:50942763:C:T:snp  10 50942763 50942763   C   T  snp OGDHL       0,0       0,0       0,0       0,0       0,0       0,0       1,0       1,0       0,0
## chr10:50943036:50943036:C:A:snp 10:50943036:50943036:C:A:snp  10 50943036 50943036   C   A  snp OGDHL       0,0       0,0       0,0       0,0       0,0       0,0       2,0       1,0       0,0
## chr10:50943058:50943058:A:T:snp 10:50943058:50943058:A:T:snp  10 50943058 50943058   A   T  snp OGDHL       0,0       0,0       0,0       0,0       0,0       0,0       3,0       2,0       0,0
##                                 03345.TAD 03345.FAD 03345.DUP
## chr10:50942493:50942493:C:T:snp       0,0       0,0       0,0
## chr10:50942558:50942558:G:A:snp       0,0       0,0       0,0
## chr10:50942763:50942763:C:T:snp       0,0       0,0       0,0
## chr10:50943036:50943036:C:A:snp       1,0       1,0       0,0
## chr10:50943058:50943058:A:T:snp       1,0       1,0       0,0

indels<-final.FAD
#the.samples<-Controls
#indels<-test.data
# alt.reads.reference.calls<-function(indels,the.samples,AD.extension="FAD",threshold=1,prefix="",suffix=""){
#   have.AD<-paste(the.samples,AD.extension,sep=".") %in% colnames(indels)
#   if(sum(have.AD)!=0){
#     allele.depths.group<-indels[,paste(the.samples,AD.extension,sep=".")]
#     genotypes<-indels[,paste(the.samples,"GT",sep=".")]
#     
#     allele.depths.group[genotypes=="1/1" | genotypes=="0/1" | genotypes=="NA" | is.na(genotypes) ] <- NA
#     
#     summary.depths.group<-apply(as.matrix(allele.depths.group),1,allele.summary.with.apply)
#     # summary.depths.group<-apply(allele.depths.group,1,allele.summary.with.threshold.apply,1)
#     summary.depths.group.threshold<-apply(as.matrix(allele.depths.group),1,allele.summary.with.threshold.apply,2)
#     
#     summary.depths.group<-t(summary.depths.group)
#     summary.depths.group.threshold<-t(summary.depths.group.threshold)
#     
#   }else{
#     print("No .AD columns in input data:alt.reads.reference.calls")
#     summary.depths.group<-matrix(data=0,nrow=dim(indels)[1],ncol=3)
#   }
#   
#   colnames(summary.depths.group)<-paste(prefix,c("ALT.reads","REF.reads","Read.Balance"),suffix,sep="")
#   rownames(summary.depths.group)<-rownames(indels)
#   
#   colnames(summary.depths.group.threshold)<-paste(prefix,c("ALT.reads","REF.reads","Read.Balance","Num.pos","Num.samples"),".thresh",suffix,sep="")
#   rownames(summary.depths.group.threshold)<-rownames(indels)
#   summary.depths.group<-cbind(summary.depths.group,summary.depths.group.threshold)
#   
#   return(summary.depths.group)
#   
# }




#indels<-final.FAD
#indels<-final.FAD[,c(1:10,69,98,99,101,310:312,1123:1128)]
#indels<-read.delim("/media/TRI-T-DRIVE-taneupan/uq/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/test.indels.txt",header=TRUE,skip=0,sep="\t",check.names=FALSE,fill=TRUE,stringsAsFactors=FALSE,strip.white=TRUE,na.strings="",quote="")
#Controls<-c("THP-1-normal-parent","10694","10695","11180" )
#cancer<- c("THP-1-AS703026resist","THP-1-Selumetinibresist")
#only use A,D columnes where D>= threshold
Controls


indels<-final.FAD
indels[1:10,1:10]
#This replaces "NA," and ",NA" to 0 
indels[] <- lapply(indels, sub, pattern = "NA(?=,[0-9])|(?<=[0-9],)NA", 
                 replacement = "0", perl = TRUE)
rownames(indels)<-indels[,"key.start"]
#write.table(indels,file="/media/TRI-T-DRIVE-taneupan/uq/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/SKDP_nextera_amanda/pvalue.indels.indels.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE, na="")
#fix(mydf)

Control.alt.counts<-alt.reads.reference.calls(indels,Controls,AD.extension="FAD",threshold=1) ## needs a sample.GT column 0/0 to start and 0/0 to end

Control.alt.counts<-alt.reads.Non.reference.calls(indels,Controls,AD.extension="FAD",threshold=1) ## needs a sample.GT column 0/1 to start and 0/1 to end
# maf.filter
# rare.in.Control 
#snp.only<-grepl("^snp",a.indel.ori[,"TYPE"]) ### the p.values codes uses "AD" so only work for SNPs
snp.only<-grepl("^snp",indels[,"TYPE"]) #
#is.flat<-grepl("flat$",indels[,"TYPE"])
#has.one.geno<-as.numeric(summary.geno.extra.ori[,"ALT.Alleles.AML"])>0
alt.count.thresh<-1  # use  p[true.sample.minor.counts < alt.count.thresh]<-1 default is less than one
#has.one.geno.ori<-has.one.geno
#to.recover<-snp.only & has.one.geno & !is.flat   ## maf.filter & rare.in.Control probably not required
to.recover<-snp.only
recover.samples<-c(cancer)
sum(to.recover)
AD.lower.tail <- FALSE # counting when more reads rather than less

#geno.p<-genotype.p.values(a.indel.stats[to.recover ,] ,c(cancer,PDs),AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh) ## COULD USE NORMAL HERE
geno.p<-genotype.p.values(indels[to.recover ,],recover.samples,AD.extension="FAD",Control.alt.counts[to.recover,"Read.Balance"]/100,cellularity,alt.count.thresh,AD.lower.tail ) ## COULD USE NORMAL HERE

# geno.p[1:10,]
# THP-1-AS703026resist THP-1-Selumetinibresist
# chr10:93616:C:T:snp:TUBB8:NM_177987:exon4:c.G716A:p.C239Y,                3.009264e-12            6.134066e-01
# chr10:93694:T:C:snp:TUBB8:NM_177987:exon4:c.A638G:p.K213R,                4.270693e-28            6.807601e-06
# chr10:1142208:T:C:snp:WDR37                                               0.000000e+00            0.000000e+00
# chr10:4703928:T:C:snp:LINC00705                                           1.000000e+00            1.000000e+00
# chr10:5202104:C:T:snp:AKR1CL1                                             0.000000e+00            0.000000e+00
# chr10:5255025:A:G:snp:AKR1C4:NM_001818:exon7:c.A749G:p.Q250R,             0.000000e+00            0.000000e+00
# chr10:5683885:A:C:snp:ASB13:NM_024701:exon5:c.T557G:p.L186R,              0.000000e+00            0.000000e+00
# chr10:5790420:T:C:snp:FAM208B:NM_017782:exon15:c.T5036C:p.V1679A,         0.000000e+00            0.000000e+00
# chr10:7601836:C:T:snp:ITIH5                                               0.000000e+00            0.000000e+00
# chr10:7763661:A:G:snp:ITIH2:NM_002216:exon8:c.A788G:p.N263S,              0.000000e+00            0.000000e+00

####################################################################################
p.threshold.z.thresh<-6
p.threshold=2*(pnorm(abs(c(p.threshold.z.thresh)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
p.threshold
found.genotype<-  geno.p <= p.threshold
geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
geno.p[!found.genotype]<-"0/0"
colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")



genotypes<-indels[to.recover ,paste(recover.samples,"GT",sep=".")] # check genotypes and geno.p in same order?
ref.call<-genotypes =="0/0"
geno.p[!ref.call]<-genotypes[!ref.call]

####################################################################################

colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")
sum(!( colnames(geno.p) %in% colnames(indels))) # must ge zero


to.transfer<-colnames(geno.p)[colnames(geno.p) %in% colnames(indels)]

posns<-match(rownames(indels),rownames(geno.p))
missing<-is.na(posns)
sum(missing)
sum(!missing)
dim(geno.p)


#######################################################
indels[!missing,to.transfer]<-geno.p[posns[!missing],to.transfer]
#######################################################
final.output<-indels[,c(colnames(indels)[1:71],wanted.samplesGT,wanted.samplesAD,wanted.samplesFAD)]

output<-cbind(final.output,geno.p)
sum(rownames(output)==final.output$key.start)  #has to be the number of rows of indels
#output<-final.output[,c(colnames(final.output)[1:71],wanted.samplesGT,wanted.samplesAD,wanted.samplesFAD)]
########### DONE

write.table(output,file="/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/pvalue_output_Aug18_0-1.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE, na="")


####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################

# 
# save(list=c("pheno.ori","a.indel.ori","Control.alt.counts","cellularity","summary.geno.extra.ori"),file="p-value.model.test.data.RData")
# 
# 
# 
# 
# 
# 
# ############# z scored based on a binomial distribution
# #README
# # http://en.wikipedia.org/wiki/Binomial_distribution
# #  Mean in n*p
# # variance np(1-p) sd =sqrt(var)
# # Models use n as the coverage or reads at a given postion
# ## Z=(X-np)/sqrt(np(1-p)) :z-score of the minor allele frequency (the number of standard deviations and observation is above the mean: (np is mean in Normal, X is measurement in tumor(minor allele counts in treatments) , the standard deviation is sqrt(np(1-p))
# # and p as the frequency of the alternative allele in the control/Normal sample.
# 
# # the programs cacluales a z-score (the number of stand deviations the allele frequnecy in the tumor sample is away from the control group)
# # the z-score is then converged to a p-value.
# # a count for the control/tumor mox ins included which extimates the number of alternative alleles might be due to the control sample along
# 
# # the aity is used to estimate the mix of tumor to normal tissue
# # alt.reads.reference.calls is used to help estimate the alternative allele frequency in control or normal samples.
# #http://www.regentsprep.org/regents/math/algtrig/ats7/blesson3.htm
# #p= normal alt reads/normal ref reads
# #n= resistant ref reads
# #x= resistant alt reads 
# #abs(c)=abs(z)
# 
p=0.001
n=100

n*p

var<- n*p*(1-p)
var

X<-4
#X<-16.5
#n<-70
#p<-0.1
Z<-(X-n*p)/sqrt(n*p*(1-p))
Z
# #z<- (true.sample.minor.counts- true.sample.cov*normal.allele.freq) / (sqrt(true.sample.cov*normal.allele.freq*(1-normal.allele.freq)))
# 
# 
2*(pnorm(Z, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
# 
# sqrt(var)
# ############################# rescue somatics with few calls with possion model
# ## bam.FADarity<-read.delim("/media/UQCCG/Sequencing/Data/Sequence_Genotypes/2015/2015-02-11_LeoPharmFinalAnalysis/Sequenza/NSLM.bam.FADarity.summary.use.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)
# 
# ## cellularity[1:5,]
# 
# 
# 
# 
# #code.dir<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/Programming/VersionControl_GitRepository/UQCCG_Pipeline_Rscripts"
# #working.directory<-("/media/UQCCG/UQCCG-Projects/Achal") 
# working.directory<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/" # directory where   "p-value.model.test.data.RData" located  
# setwd(working.directory)
# source("annotate_SNPs_subroutines.r")
# options("width"=250,"max.print"=1000)
# core.ann<-c("chr","start","end","REF","ALT","TYPE")
# 
# #working.directory<-"/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal" # directory where   "p-value.model.test.data.RData" located  
# #setwd(working.directory)
# 
# #project.file<-"2015-03-25_SKDPNexteraRun_CSV.THP.wanted.All-maf-filtered_PAul.csv"
# project.file<-"test2015-03-25_SKDPNexteraRun.THP.wanted.All-maf-filtered.csv"
# 
# cellularity<-read.delim("cellularity.csv",header=T,sep="\t",fill=TRUE,skip=0,stringsAsFactors=FALSE)
# cellularity[1:5,]
# 
# column.labels<-read.delim(project.file,header=F,nrows=1,sep="\t",fill=TRUE,stringsAsFactors=FALSE,na.strings="",quote="\"")
# num.vars<-dim(column.labels)[2]
# a.indel<-scan(project.file,what=character(num.vars),skip=1,sep="\t",fill=TRUE,na.strings="",quote="\"")
# num.lines<-length(a.indel)/(num.vars)
# dim(a.indel)<-c(num.vars,num.lines)
# a.indel<-t(a.indel)
# colnames(a.indel)<-column.labels
# ########################################
# the.chr<-unique(a.indel[,"chr"])
# print(paste("Doing Chromosome ",the.chr))
# 
# if(sum(!grepl("^chr",the.chr))>0){
#   a.indel[,"chr"]<-paste("chr",a.indel[,"chr"],sep="")
# }
# ###################################### a.indel!!
# a.indel[1:5,"chr"]
# key<-build.key(a.indel,core.ann)
# summary.geno.extra<-a.indel
# rownames(summary.geno.extra)<-key
# rownames(a.indel)<-key
# 
# #load("p-value.model.test.data.RData")
# 
# #ls() ## see what you loaded
# # a.indel<-a.indel.ori
# # summary.geno.extra<-summary.geno.extra.ori
# 
# #normals<-pheno.ori[pheno.ori[,"normal"],"SAMPLE"]
# #normals<-normals[normals!="LPH-001-27_PD"]
# #normals
# test<-a.indel
# #colnames(test)[90:92]<-c("THP-1-AS703026resist.ReCal.sort.bam.GT","THP-1-normal-parent.ReCal.sort.bam.GT","THP-1-Selumetinibresist.ReCal.sort.bam.GT")
# a.indel <-cbind(a.indel, "THP-1-normal-parent.ReCal.sort.bam.GT"=a.indel[,"THP-1-normal-parent.GT"], "THP-1-AS703026resist.ReCal.sort.bam.GT"= a.indel[,"THP-1-AS703026resist.GT"],"THP-1-Selumetinibresist.ReCal.sort.bam.GT"= a.indel[,"THP-1-Selumetinibresist.GT"])
# 
# 
# #Control<-c("THP-1-normal-parent")
# #Control
# Control<-c("THP-1-normal-parent.ReCal.sort.bam")
# #Control.alt.counts <- alt.reads.reference.calls(a.indel,Control,AD.extension="bam.AD",threshold=1) 
# Control.alt.counts<-alt.reads.Non.reference.calls(test,Control,AD.extension="bam.AD",threshold=1) # 0/1 start
# 
# treatments<-c("THP-1-AS703026resist.ReCal.sort.bam.GT","THP-1-Selumetinibresist.ReCal.sort.bam.GT")
# #treatments<-c("THP-1-AS703026resist","THP-1-Selumetinibresist")
# treatments
# 
# 
# ########################################
# ########################################
# ########################################
# 0/0 <->  0/0 P value different
# 0/1 <->  0/1 
# 
# 
# 
# snp.only<-grepl("^snp",a.indel[,"TYPE"]) ### the p.values codes uses "AD" so does ONLY work for SNPs
# 
# #has.one.geno<-as.numeric(a.indel[,"ALT.Alleles.cancer"])>0
# 
# geno.p<-genotype.p.values(a.indel[snp.only,],treatments,Control.alt.counts[snp.only, "Read.Balance"]/100,cellularity,AD.extension="bam.AD")
# geno.p
# 
# #write.table(geno.p,file="/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/geno.p_0_0_nonref_0_0.txt",col.names=TRUE,row.names=TRUE,sep="\t",quote=FALSE, na="")
# 
# #p.threshold=0.0026 # z=3:   2*(pnorm(abs(c(1:8)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
# p.threshold=2*(pnorm(abs(c(6)), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE))
# p.threshold
# 
# #geno.p[geno.p>p.threshold]<-"mor.thold"
# 
# geno.p <- cbind(geno.p,Control.alt.counts[snp.only, ])
# geno.p[1:10,]
# 
# snp.only.indel <- (a.indel[snp.only,])
# 
# final <- cbind (snp.only.indel[,1:92], geno.p, snp.only.indel[,93:dim(snp.only.indel)[2]])
# 
# #final<-final[!grepl("flat$",final[,"TYPE"]),]
# 
# 
# write.table(final,file="/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/0_1_ref_0_1_with_flat.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE, na="")
# 
# write.table(test,file="/media/TRI-T-DRIVE-taneupan/uqdi/Core_Services/UQCCG/UQCCG-Projects/Achal/Pvalue_enrichment_july9/0_1_ref_0_1_with_flat.txt",col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE, na="")
# 
# 
# found.genotype <- geno.p<=p.threshold
# #geno.p[found.genotype]<-"0/1" # p.threshold<- 0.0026
# #geno.p[!found.genotype]<-"0/0"
# #colnames(geno.p)<-paste(colnames(geno.p),".GT",sep="")
# 
# 
# ########################################################################################################################
# 
# x <- as.matrix(read.csv(text="start,A,T,G,C,REF,ALT,TYPE 
# chr20:5363934,95,29,14,59,C,T,snp
# chr5:8529759,,,,,G,C,snp
# chr14:9620689,65,49,41,96,T,G,snp
# chr18:547375,94,1,51,67,G,C,snp
# chr8:5952145,27,80,25,96,T,T,snp
# chr14:8694382,68,94,26,30,A,A,snp
# chr16:2530921,49,15,79,72,A,T,snp:2530921
# chr16:2530921,49,15,79,72,A,G,snp:2530921
# chr16:2530921,49,15,79,72,A,T,snp:2530921flat
# chr16:2530331,9,2,,,A,T,snp:2530331
# chr16:2530331,9,2,,,A,G,snp:2530331
# chr16:2530331,9,2,,,A,T,snp:2530331flat
# chr16:2533924,42,13,19,52,G,T,snp:flat
# chr16:2543344,4,13,13,42,G,T,snp:2543344flat
# chr16:2543344,42,23,13,42,G,A,snp:2543344
# chr14:4214117,73,49,18,77,G,A,snp
# chr4:7799768,36,28,1,16,C,A,snp
# chr3:9141263,27,41,93,90,A,A,snp", stringsAsFactors=FALSE))
# 
# 
# ##################################
# 
# df <- as.data.frame(x)
# types <- c("A", "T", "G", "C")
# df[types][is.na(df[types])] <- 0
# head(newBamAD(df))
# 
# newBamAD <- function (x,base.types=c("A","C","G","T")) {
# 
#   # the version above
#   rownames(x) <- 1:nrow(x)
#   ref <- x[cbind(1:nrow(x), x[, "REF"])]
#   alt <- x[cbind(1:nrow(x), x[, "ALT"])]
#   which.flat <- grep("flat$", x[, "TYPE"])
#   
#   alt[which.flat] <-  sapply(which.flat, function (i,base.types) {
#     sum(as.numeric(x[i, c( base.types[!( base.types %in% x[i, "REF"])] )] ) ,na.rm=TRUE) },base.types)
#   cbind(x[,c("start","REF","ALT","TYPE")],bam.AD=paste(ref, alt, sep=','))
#   # cbind(x, bam.AD=paste(ref, alt, sep=','))
# }
# newBamAD(df)
