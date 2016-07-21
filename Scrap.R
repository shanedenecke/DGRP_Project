library(ggplot2)
library(plyr)

File_Name=file.choose()
mainDir=dirname(File_Name)
setwd(mainDir)
Data=read.csv(File_Name)
Data=Data[complete.cases(Data),]

Corrected_Data=data.frame()
for (a in as.character(unique(Data$Rep))){; for (b in as.character(unique(Data$Metabolite))){; for(c in as.character(unique(Data$Genotype))){
  sub=subset(Data,Rep==a & Metabolite==b & Genotype==c)
  if(length(sub$Rep)>0){
		sub$Ratio_to_Media=sub$Larval_Amount/sub$Media_Amount
		Corrected_Data=rbind(Corrected_Data,sub)
	}else{
		print("Such is Life")
		}}
}}

CI=function(x){
  y=qnorm(.975)*(sd(x)/sqrt(length(x)))
  return(y)
}

Summary=ddply(Corrected_Data, c("Metabolite","Group","Rank","Genotype","X25ppm_WI","X100ppm_WI"), summarise,
                   N    = length(Larval_Amount),
                  Mean_Larval_Amount = mean(Larval_Amount),
                  Sd_Larval_Amount   = sd(Larval_Amount),
                  Se_Larval_Amount   = Sd_Larval_Amount/sqrt(N),
                  CI_Larval_Amount   = qnorm(0.975)*Se_Larval_Amount,
                  Mean_Media_Amount = mean(Media_Amount),
                  Sd_Media_Amount   = sd(Media_Amount),
                  Se_Media_Amount   = Sd_Media_Amount/sqrt(N),
                  CI_Media_Amount   = qnorm(0.975)*Se_Media_Amount,
                  Mean_Ratio = mean(Ratio_to_Media),
                  Sd_Ratio   = sd(Ratio_to_Media),
                  Se_Ratio   = Sd_Ratio/sqrt(N),
                  CI_Ratio   = qnorm(0.975)*Se_Ratio)
Summary$Sd_Larval_Amount=NULL
Summary$Sd_Media_Amount=NULL
Summary$Sd_Ratio=NULL
Summary$Se_Larval_Amount=NULL
Summary$Se_Media_Amount=NULL
Summary$Se_Ratio=NULL
             
             

attach(Summary)
for (M in unique(Metabolite)){
  subT=subset(Summary,Metabolite==M)
  p=ggplot(data=subT,aes(Genotype,c(Mean_Larval_Amount,Mean_Media_Amount),fill=)) 
  p=p+geom_bar(stat="identity",position="dodge")
  p=p+geom_errorbar(stat="identity",color="black")
  
  
  
  	x=ggplot(data=GLM_Summary, aes(x=Genotype, y=Mean,fill=Genotype)) +
			geom_errorbar(aes(ymin=Mean,ymax=Mean+CI), width=.5,size=.5,position=position_dodge()) +
			geom_bar(stat="identity") +
			geom_bar(stat="identity",color="black", show_guide=FALSE) +
			xlab("Genotypes") +
			ylab(expression(bold(paste("-",beta," values",sep="")))) +
			scale_fill_manual(values=c(Color_Base_Genotype)) +
			theme(text=element_text(face="bold",family="Times")) +
			theme(axis.text.x=element_text(size=8,face="bold",family="Times")) +
			theme(text=element_text(size=10,face="bold",family="Times"),plot.title=element_text(size=16,face="bold",family="Times"),legend.title=element_text(size=12,face="bold",family="Times"))+
			background_grid(major = "xy", minor = "none",colour.major="grey75")+	
			ylim(0,lim) +
			ggtitle(paste(D," ",I,sep=""))
			save_plot(filename=paste(D,I,"GLM.pdf",sep=" "),plot=x)	
}}




# Get Larvae/Media Ratios from respective values
Corrected_Data=data.frame()
for (a in as.character(unique(Data$Rep))){; for (b in as.character(unique(Data$Metabolite))){; for(c in as.character(unique(Data$Genotype))){
  sub=subset(Data,Rep==a & Metabolite==b & Genotype==c)
  if(length(sub$Rep)>3){
		sub.sub <- sub[1,1:5]
        sub.sub$Measurement_Type="Ratio"
        sub.sub$Measurement_Value=sub$Measurement_Value[which(sub$Measurement_Type=="Larvae")]/sub$Measurement_Value[which(sub$Measurement_Type=="Media")]
        Corrected_Data=rbind(Corrected_Data,sub.sub)
        #sub$Ratio=sub$Larvae/sub$Media
		#Corrected_Data=rbind(Corrected_Data,sub)
	}else{
		print("Such is Life")
		}}
}}

Corrected_Data=rbind(Corrected_Data,Data)

## In Base Plots


for (m in Metabolites){
  for (p in Comparison){
    pdf(file=paste(m,p,"dotplot.pdf",sep="_"))
    sub=subset(Corrected_Data,Metabolite==m)
    submat=as.matrix(sub)
    par(mar=c(2,2,2,2),mgp=c(3,1,0),mfrow=c(1,1),xpd=TRUE)
    dotchart(as.numeric(submat[,p]),group=sub$Genotype,gcolor="black",lcolor=NA,color=dot.col)
  }}


## Rects stuff
back.data <- data.frame(xstart=c(0,9),xend=c(9,18),col=c("blue","red"))
    
    rects <- data.frame(xstart = seq(0,80,20), xend = seq(20,100,20), col = letters[1:5])

ggplot() + 
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf, fill = col), alpha = 0.4) +
  geom_line(data = dat, aes(x,y))


Make Boxplots for your Single Metabolite Stuff



    
    
    
    
   
    
    
    
    
     gp=gp+geom_rect(data=sub, aes(xmin=0, xmax=10,ymin=-Inf,ymax=Inf,fill="red", alpha = 0.4))
    
    
    
pdf(file=paste(m,p,"boxplot.pdf",sep="_"))
par(mar=c(7,6,4,2),mgp=c(3,1,0),mfrow=c(1,1),xpd=TRUE)
boxplot(formula=get(p)~Genotype,data=Corrected_Data,subset=Metabolite==m,col=total.col,axes=F,main=paste(m,p,sep=" "),
staplecol=total.col,whiskcol=total.col,staplelwd=1.4,whisklwd=1.4,names=c("1,2,3")) 
axis(side=1,at=1:20,labels=F,cex=1.2,family="serif",tick=T)
axis(side=2,at=seq(0,signif(max(subset(Corrected_Data,Metabolite==m)[,p]),2),length.out=6),font=2,cex=1.2,las=1,family="serif")
mtext(side=1, text="Susceptible",at=5,font=2,cex=1.5,line=5,family="serif")
mtext(side=1, text="Resistant",at=15,font=2,cex=1.5,line=5,family="serif")
mtext(side=2, text=paste("Concentration of ",m," (ppb)",sep=""),family="serif",font=2,cex=1.5,line=3.5)
text(seq(par("usr")[1]+1.5,par("usr")[2]-1,length.out=20), par("usr")[3]-0.5, srt = 45, xpd = TRUE,adj=1.3,labels =as.character(unique(Corrected_Data$Genotype)), font=2, family="serif",cex=.8)
box(bty="l")
dev.off()
}}

## Summarize your Data for other plots
Summary=ddply(Corrected_Data, c("Metabolite","Group","Rank","Genotype","X25ppm_WI","X100ppm_WI"), summarise,
                   N    = length(Larvae),
                  Mean_Larvae = mean(Larvae),
                  Se_Larvae   = sd(Larvae)/sqrt(N),
                  CI_Larvae   = qnorm(0.975)*Se_Larvae,
                  Mean_Media = mean(Media),
                  Se_Media   = sd(Media)/sqrt(N),
                  CI_Media   = qnorm(0.975)*Se_Media,
                  Mean_Ratio = mean(Ratio),
                  Se_Ratio   = sd(Ratio)/sqrt(N),
                  CI_Ratio   = qnorm(0.975)*Se_Ratio)

Summary=arrange(Summary,Rank)



## Single Metabolite BarPlots. Axis still giving me problems

for (m in Metabolites){
  for (p in Comparison){
    sub.m=subset(Summary,Metabolite==m)
    sub.graph=sub.m[,c("Genotype","Group","Rank","X25ppm_WI","X100ppm_WI",paste("Mean_",p,sep=""),paste("Se_",p,sep=""),paste("CI_",p,sep=""))]
    
    attach(sub.graph)
    mean=get(paste("Mean",p,sep="_"))
    se=get(paste("Se",p,sep="_")) 

pdf(file=paste(m,p,"barplot.pdf",sep="_"),width=10)
par(mar=c(6,6,3,3),mgp=c(3.2,1,0),mfrow=c(1,1))
bp=barplot(mean,col=total.col,ylim=c(0,3+max(mean+se)),family="serif",
ylab=list("Metabolite Level",font=2,cex=1.8),main=m,cex.main=2,axes=F,las=2)
segments(bp,mean, bp,mean+se,col="black",lwd=2)
segments(x0=bp-.3, y0=mean+se,x1=bp+.3,lwd=2,col="black")
axis(side=2,at=seq(0,signif(max(subset(Corrected_Data,Metabolite==m)[,p]),2),length.out=6),font=2,cex=1.2,las=1,family="serif")
text(seq(par("usr")[1]+1.5,par("usr")[2]-1.2,length.out=20), par("usr")[3]-0.25, srt = 45, xpd = TRUE,adj=1.2,labels =paste(unique(Genotype)), font=2, family="serif",cex=.8)
mtext(side=1, text="Susceptible",at=5,font=2,cex=1.5,line=3.5,family="serif")
mtext(side=1, text="Resistant",at=15,font=2,cex=1.5,line=3.5,family="serif")
dev.off()
  }
}

## Group WI value at 25ppm with individual metabolites

for (m in Metabolites){
  for (p in Comparison){
    sub.m=subset(Summary,Metabolite==m)
    sub.graph=sub.m[,c("Genotype","Group","Rank","X25ppm_WI","X100ppm_WI",paste("Mean_",p,sep=""),paste("Se_",p,sep=""),paste("CI_",p,sep=""))]
    sub=cbind(sub.graph[,paste("Mean_",p,sep="")],sub.graph$X25ppm_WI)
    rownames(sub)=as.character(unique(sub.graph$Genotype))
    colnames(sub)=c(p,"WI_25")
    sub=t(sub)
    par(mar=c(6,6,3,3),mgp=c(3.2,1,0),mfrow=c(1,1))
    barplot(sub,beside=T)



a=c(



attach(Summary)
for (M in unique(Metabolite)){
  subT=subset(Summary,Metabolite==M)
  p=ggplot(data=subT,aes(Genotype,c(Mean_Larval_Amount,Mean_Media_Amount),fill=)) 
  p=p+geom_bar(stat="identity",position="dodge")
  p=p+geom_errorbar(stat="identity",color="black")
  
  
  
  	x=ggplot(data=GLM_Summary, aes(x=Genotype, y=Mean,fill=Genotype)) +
			geom_errorbar(aes(ymin=Mean,ymax=Mean+CI), width=.5,size=.5,position=position_dodge()) +
			geom_bar(stat="identity") +
			geom_bar(stat="identity",color="black", show_guide=FALSE) +
			xlab("Genotypes") +
			ylab(expression(bold(paste("-",beta," values",sep="")))) +
			scale_fill_manual(values=c(Color_Base_Genotype)) +
			theme(text=element_text(face="bold",family="Times")) +
			theme(axis.text.x=element_text(size=8,face="bold",family="Times")) +
			theme(text=element_text(size=10,face="bold",family="Times"),plot.title=element_text(size=16,face="bold",family="Times"),legend.title=element_text(size=12,face="bold",family="Times"))+
			background_grid(major = "xy", minor = "none",colour.major="grey75")+	
			ylim(0,lim) +
			ggtitle(paste(D," ",I,sep=""))
			save_plot(filename=paste(D,I,"GLM.pdf",sep=" "),plot=x)	
		}}
