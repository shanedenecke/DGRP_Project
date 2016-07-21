## Import Libraries
library(qqman)
library(plyr)
manhattan2 <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
    "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
    genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
    ...) 
{
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x))) 
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]])) 
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]])) 
        d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }
    else {
        d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        options(scipen = 999)
        d$pos = d$BP/1e+06
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position(Mb)")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index == 
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
                  lastbase
            }
            ticks = c(ticks, (min(d[d$CHR == i, ]$pos) + max(d[d$CHR == 
                i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
        las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
            ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
        names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos, 
                logp, col = col[icol], pch = 20, ...))
            icol = icol + 1
        }
    }
    if (suggestiveline) 
        abline(h = suggestiveline, col = "lightblue")
    if (genomewideline) 
        abline(h = genomewideline, col = "darkblue")
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "red", pch = 20, 
            ...))
    }
}


## Set Working Directory and Import Data
Input.Dir=dirname(file.choose())
Output.Dir=dirname(file.choose())
Doses=c("25ppm","100ppm")


## Run your ultimate for loop that will do everything for 2 doses
for (d in Doses){
setwd(Input.Dir)
## Apply bonferroni and FDR corrections to Association Data
Initial=read.table(paste("gwas.all_Mean_",d,".assoc",sep=""),nrow=100)
classes=sapply(Initial,class)
Data=read.table(paste("gwas.all_Mean_",d,".assoc",sep=""),sep=" ",header=T,comment.char="",colClasses=classes)
Data$SNP=as.character(Data$ID)
q=p.adjust(Data$SingleMixedPval,method="fdr")
b=p.adjust(Data$SingleMixedPval,method="bonferroni")
Data$q.value=q
Data$bonferroni=b
Data=arrange(Data,SingleMixedPval)
write.csv(Data,file=paste("pAdjust_",d,".csv",sep=""))

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
png(file=paste("Manhattan_Plot_",d,".png",sep=""),width=1000,height=1000)
par(mar=c(8,8,8,8),mgp=c(4,1,0),mfrow=c(1,1),family="serif",font=2)
manhattan2(Data,col=c("black","grey"),genomewideline=-log10(.05/length(Data$SNP)), 
highlight=SNP.Significant,main=d,ylim=c(0,-log10(1e-8)),chrlabs=c('X','2L','2R','3L','3R'),
cex.main=4,cex.axis=2,cex.lab=3.5,font.lab=2)
#mtext(side=2,text=expression(bold(paste("-log"[10]*"(p)"))),font=2,cex=3,line=3,family="serif")
#mtext(side=1,text="Chromosome",font=2,cex=3,line=3,family="serif")
#axis(side=1,at=1:5,labels=,tick=TRUE,cex.axis=2.5,font=2)
#axis(side=2,at=1:7,labels=TRUE,tick=TRUE,cex.axis=2.5,padj=.5)
box(bty="l")
dev.off()   

png(paste(file="QQplot_",d,".png",sep=""),width=800,height=800) 
par(mar=c(7,7,7,7),mgp=c(2,1.5,0),mfrow=c(1,1),family="serif",font=2)
qq(Data$P, main = paste("Q-Q plot of ",d,sep=""),
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
