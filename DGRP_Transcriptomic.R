## Import Data and Set basic parameters
library(qqman)
library(dplyr)
library(data.table)
setwd("/home/sdenecke/Dropbox/PhD_Projects/DGRP_Paper/Results/Post_GWAS/Transcriptome_Associations")
raw.phen.data <- fread("Schmidt_DDT_Phenotype.csv",header=T,sep=",")
f.trans.data <- fread("F185TRANS_ANNOT.csv", header=T,sep=",")
m.trans.data <- fread("M185TRANS_ANNOT.csv", header=T,sep=",")
comb.trans.data <- rbindlist(list(m.trans.data,f.trans.data),use.names=TRUE)
Sex <- c("m","f")
Doses <- unique(raw.phen.data$Dose)


dir.create(paste("./",format(Sys.Date(),"%d.%b.%y")))
setwd(paste("./",format(Sys.Date(),"%d.%b.%y")))

 
    ## Get the genotypes for your transcriptome data and your phenotype data  
    trans.gen <- colnames(f.trans.data[-1:-3])
    phen.gen <- as.character(raw.phen.data$Genotype)
    common.gen <- intersect(trans.gen,phen.gen)

    phen.data <- raw.phen.data[phen.gen %in% common.gen,]
    trans.data <- cbind(comb.trans.data[,1:3,with=F],
                        comb.trans.data[,trans.gen %in% common.gen,with=F])
    trans.data$sex <- with(comb.trans.data,c(rep("m",length(FBID)/2),
                                                  rep("f",length(FBID)/2))) 
  
for (d in Doses){
      for (s in Sex){
            sub.trans <- subset(trans.data,sex==s)
            symbol <- sub.trans$SYMBOL
            sub.trans[,c("FBID","SYMBOL","LOCATION","sex"):=NULL]
            test.trans <- sapply(select(sub.trans,starts_with("RAL")),as.numeric)
            sub.phen <- as.numeric(subset(phen.data,Dose==d)$Phenotype)
      # Fit a linear model to each transcript
    output=matrix(NA,nrow=length(test.trans[,1]),ncol=3)
    colnames(output)=c("estimate","pval","adjrsq")
    for (i in 1:nrow(test.trans)) {
	  test <- summary(lm(sub.phen~test.trans[i,]))
	  est <- test$coefficients[2,1]
	  pval <- test$coefficients[2,4]
	  rsq <- test$adj.r.squared
  	  output[i,] <- c(est,pval,rsq)
    }
      

    ## Clean up and write your data to a .csv file
    rownames(output)=sub.trans[,"SYMBOL"]
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



















