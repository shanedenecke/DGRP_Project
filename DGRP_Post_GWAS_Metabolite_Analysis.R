
## Import Libraries and Data into R
library(plyr)
library(ggplot2)
File_Name=file.choose()
mainDir=dirname(File_Name)
setwd(mainDir)
Data=read.csv(File_Name)
Data=Data[complete.cases(Data),]

##Create a function that calculates 95% confidence intervals

CI=function(x){
  y=qnorm(.975)*(sd(x)/sqrt(length(x)))
  return(y)
}


# Get Larvae/Media Ratios from respective values
Corrected_Data=data.frame()
for (a in as.character(unique(Data$Rep))){; for (b in as.character(unique(Data$Metabolite))){; for(c in as.character(unique(Data$Genotype))){
  sub=subset(Data,Rep==a & Metabolite==b & Genotype==c)
  if(length(sub$Rep)>0){
		sub$Ratio=sub$Larvae/sub$Media
		Corrected_Data=rbind(Corrected_Data,sub)
	}else{
		print("Such is Life")
		}}
}}

Corrected_Data=arrange(Corrected_Data,Rank)
Corrected_Data$Genotype=factor(Corrected_Data$Genotype,levels=unique(Corrected_Data$Genotype))
Metabolites=unique(Corrected_Data$Metabolite)
Comparison=c("Media","Larvae")

## Get Colour pallete for resistance and susceptible genotypes
res.col=colorRampPalette(c("orangered1","orangered4"))(9)
sus.col=colorRampPalette(c("steelblue1","steelblue4"))(9)
total.col=c(sus.col,res.col)
dot.col=c()
for(a in total.col){
  if(paste(unique(Corrected_Data$Genotype)[which(total.col==a)])=="RAL409"){
    dot.col=c(dot.col,rep(a,3))
  }else{
  dot.col=c(dot.col,rep(a,4))
  }
}

## In Base Plots


for (m in Metabolites){
  for (p in Comparison){
    pdf(file=paste(m,p,"dotplot.pdf",sep="_"))
    sub=subset(Corrected_Data,Metabolite==m)
    submat=as.matrix(sub)
    par(mar=c(2,2,2,2),mgp=c(3,1,0),mfrow=c(1,1),xpd=TRUE)
    dotchart(as.numeric(submat[,p]),group=sub$Genotype,gcolor="black",lcolor=NA,color=dot.col)
  }}

back.data <- data.frame(xstart=c(0,9),xend=c(9,18),col=c("blue","red"))
    
    rects <- data.frame(xstart = seq(0,80,20), xend = seq(20,100,20), col = letters[1:5])


#As Baptiste points out, the order of the geom's matters, so putting your data as last will 
#make sure that it is plotted "on top" of the background rectangles. Updated code, but
#did not update the JPEG...I think you'll get the point.

ggplot() + 
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), alpha = 0.4) +
  geom_line(data = dat, aes(x,y))


for (m in Metabolites){
  for (p in Comparison){
    sub=subset(Corrected_Data)
    pdf(file="test.pdf")
    gp=ggplot(sub,aes_string(x='Genotype',y=p)) 
   # gp=gp+geom_rect(aes(xmin=c(0,9),xmax=c(9,18),ymin=-Inf, ymax = Inf, fill = "red","blue"), alpha = 0.4)
    gp=gp+geom_dotplot(binaxis='y', stackdir='center',binwidth=1,colour="black",dotsize=2,fill=dot.col) 
    gp=gp+ylab(paste("Concentration of ",m," [ppb]",sep=""))
    gp=gp+stat_summary(geom="errorbar",colour='black',size=1.5,width=.6)
    gp=gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.6,colour='black',size=1.25)
    gp=gp+ggtitle(paste(m,p))
    gp=gp+theme(text=element_text(size=18,family="Helvetica-Narrow",face="bold"),
      axis.text.x=element_text(angle = 45, hjust = 1,size=16,colour="black"),
      axis.text.y=element_text(size=14,family="",face="bold"),
      panel.background=element_rect(fill="grey99"),
      panel.border=element_rect(colour="black",fill=NA),
      panel.grid.minor=element_blank(),panel.grid.major.x=element_blank())
    gp
    #gp=gp+geom_rect(xmin=0, xmax=10,ymin=-Inf,ymax=Inf,fill="blue",alpha)
    dev.off()
  }}


