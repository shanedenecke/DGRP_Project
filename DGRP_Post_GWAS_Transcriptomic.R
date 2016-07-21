## Import Data and Set basic parameters
library(qqman)
File_Name <- file.choose()
mainDir <- dirname(File_Name)
setwd(mainDir)
raw.data <- read.csv(File_Name)
sex <- c("M","F")
Doses <- c("25ppm","100ppm")
dir.create(file.path(getwd(),Sys.Date()))
setwd(file.path(getwd(),Sys.Date()))

## Run Two nested for loops that will give you sex based and dose based outputs

for (d in Doses){ 
  for (s in sex){

    ## Import your transcriptome data and phenotypic data
    raw.phen.data=subset(raw.data,Dose==d)
    if (s=="F") {
	    raw.trans.data=read.csv("F185TRANS_ANNOT.csv", header=TRUE)
    } else {
  	  raw.trans.data=read.csv("M185TRANS_ANNOT.csv", header=TRUE)
    }
 
    ## Get the genotypes for your transcriptome data and your phenotype data  
    trans.gen=colnames(raw.trans.data[-1:-3])
    phen.gen <- as.character(raw.phen.data$Genotype)
    common.gen<- unique(c(trans.gen,phen.gen))

    ## Return genotypes common to both datasets
    v <- character()
    for (a in common.gen){
        if (a %in% trans.gen & a %in% phen.gen){
            v <- c(v,a)
        }
    }
    phen.data <- raw.phen.data[phen.gen %in% v,]
    trans.data <- cbind(raw.trans.data[1:3],raw.trans.data[-1:-3][,trans.gen %in% v])
    
  
    # Fit a linear model to each transcript
    output=matrix(NA,nrow=0,ncol=3)
    colnames(output)=c("estimate","pval","adjrsq")
    for (i in 1:nrow(trans.data)) {
	  test=lm(as.numeric(phen.data[,"Mean"])~as.numeric(trans.data[-1:-3][i,]))
	  est=summary(test)$coefficients[2,1]
	  pval=summary(test)$coefficients[2,4]
	  rsq=summary(test)$adj.r.squared
	  out.line=c(est,pval,rsq)
  	  output=rbind(output,out.line)
    }


    ## Clean up and write your data to a .csv file
    rownames(output)=trans.data[,"SYMBOL"]
    summ=output[order(output[,"pval"]),]
    fdr=p.adjust(summ[,"pval"],method="fdr")
    bonferroni=p.adjust(summ[,"pval"],method="bonferroni")
    summ <- cbind(summ,fdr)
    summ <- cbind(summ,bonferroni)
    write.csv(summ,file=paste("TWAS_output_",s,d,".csv",sep=""))
    
    ## Make a qqplot of your transcriptome data
    png(paste(file="QQplot_",d,s,".png",sep=""),width=800,height=800) 
    par(mar=c(7,7,7,7),mgp=c(2,1.5,0),mfrow=c(1,1),family="serif",font=2)
    qq(summ[,"pval"], main = paste("Q-Q plot of ",d,s,sep=""),
    cex.main=2,cex.axis=3,axes=F,
    ylim=range(0,7),xlim=range(0,7),
    xlab="",ylab="")
    mtext(side=1,text=expression(bold(paste("Expected -log"[10]*"(p)"))),cex=3,line=5)
    mtext(side=2,text=expression(bold(paste("Observed -log"[10]*"(p)"))),cex=3,line=3)
    axis(side=1,at=1:7,labels=TRUE,tick=TRUE,cex.axis=2.5,font=2)
    axis(side=2,at=1:7,labels=TRUE,tick=TRUE,cex.axis=2.5,font=2,padj=.5)
    box(bty="l")
    dev.off()  
  } 
}




















