## Import Libraries
library(qqman)
library(plyr)

## Set Working Directory and Import Data
Input.Dir=dirname(file.choose())
Output.Dir=dirname(file.choose())
Doses=c("25ppm","100ppm")
Models=c("Mean")#,"BLUP")


## Run your ultimate for loop that will do everything for 2 doses
for (d in Doses){ 
  for (m in Models){
setwd(Input.Dir)
## Apply bonferroni and FDR corrections to Association Data
Data=read.table(paste("gwas.all_",m,"_",d,".assoc",sep=""),sep=" ",header=T)
Data$SNP=as.character(Data$ID)
q=p.adjust(Data$SingleMixedPval,method="fdr")
b=p.adjust(Data$SingleMixedPval,method="bonferroni")
Data$q.value=q
Data$bonferroni=b
Data=arrange(Data,SingleMixedPval)
write.csv(Data,file=paste("pAdjust_",m,"_",d,".csv",sep=""))

## Get Rid of all that crap and lets just leave you with what you need for the manhattan plots

## Only Run if you are testing
#Data=read.csv("pAdjust.csv")
#Data$SNP=as.character(Data$SNP)
## Only Run if you are testing



Data$MinorAllele=NULL  
Data$MajorAllele=NULL  
Data$MinorAlleleCount=NULL  
Data$MajorAlleleCount=NULL  
Data$SinglePval=NULL  
Data$ID=NULL  
Data$MAF=NULL  
Data$q.value=NULL  
Data$bonferroni=NULL  

 
## Out of those messy SNP labels we'll get you a cleaner version of your chromosme and base annotations

SNP2=unlist(strsplit(Data$SNP,split="_"))
Chrome=SNP2[seq(1, length(SNP2), 3)]
Base=as.numeric(SNP2[seq(2, length(SNP2), 3)])
Chrome=gsub('X',1,Chrome)
Chrome=gsub('2L',2,Chrome)
Chrome=gsub('2R',3,Chrome)
Chrome=gsub('3L',4,Chrome)
Chrome=gsub('3R',5,Chrome)
Chrome=as.numeric(Chrome)
Data$CHR=Chrome
Data$BP=Base
colnames(Data)=c('P','SNP','CHR','BP')

SNP.Significant=Data$SNP[which(Data$P<=.05/length(Data$SNP))]
setwd(Output.Dir)

# Ok now lets generate some mahattan plots
Data=Data[complete.cases(Data),]

png(file=paste("Manhattan_Plot_",m,"_",d,".png",sep=""),width=800,height=800)
par(mar=c(5,6,2,1),mgp=c(3,1,0),mfrow=c(1,1))
manhattan(Data,col=c("black","orange"),genomewideline=-log10(.05/length(Data$SNP)),chrlabs=c("X",'2L','2R','3L','3R'),
highlight=SNP.Significant,main=paste(m,d,sep=" "),cex.main=2,ylim=c(0,-log10(1e-8)),family="serif",font=2,cex.lab=1.5,cex.axis=1.5,
xlab="",ylab=list(expression(bold(paste("-log"[10]*"(p)"))),font=2))  
dev.off()    

png(paste(file="QQplot_",m,"_",d,".png",sep=""),width=800,height=800) 
par(mar=c(5,6,2,1),mgp=c(3,1,0),mfrow=c(1,1),bty="l") 
qq(Data$P, main=paste(m,d,sep=" "), family="serif",cex.main=2,font.axis=2,cex.lab=1.5,cex.axis=1.5,
ylab=list(expression(bold(paste("Observed -log"[10]*"(p)"))),font=2), 
xlab=list(expression(bold(paste("Expected -log"[10]*"(p)"))),font=2)) 
dev.off()

  }
}




ylab=list(expression(bold(paste("Observed -log"[10]*"(p)"))),font=2)
xlab=list(expression(bold(paste("Expected -log"[10]*"(p)"))),font=2)


