
## Import Libraries and Data into R
library(plyr)
library(ggplot2)
File_Name=file.choose()
mainDir=dirname(File_Name)
setwd(mainDir)
Data=read.csv(File_Name)
Data=Data[complete.cases(Data),]
dir.create(paste("./",format(Sys.Date(),format="%d %b %y"),sep=""))
setwd(paste("./",format(Sys.Date(),format="%d %b %y"),sep=""))


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
Comparison=c("Media","Larvae","Ratio")

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

for (m in Metabolites){
  for (p in Comparison){
    sub=subset(Corrected_Data,Metabolite==m)
    lim <- ceiling(max(sub[["Larvae"]],sub$Media))
    pdf(file=paste(m,p,"pdf",sep="."))
    gp=ggplot(sub,aes_string(x='Genotype',y=p)) 
    gp=gp+geom_dotplot(binaxis='y', stackdir='center',colour="black",dotsize=.5,fill=dot.col) 
    gp=gp+ylab(paste("Concentration of ",m," [ppb]",sep=""))
    gp=gp+stat_summary(geom="errorbar",colour='black',size=1,width=.6,fun.data="mean_se")
    gp=gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.6,colour='black',size=.8)
    gp=gp+ggtitle(paste(m,p))
    gp <- gp+ylim(c(0,lim))
    gp=gp+theme(text=element_text(size=18,family="",face="bold"),
      axis.text.x=element_text(angle = 45, hjust = 1,size=16,colour="black"),
      axis.text.y=element_text(size=14,family="",face="bold"),
      panel.background=element_rect(fill=NA),
      panel.border=element_rect(colour="black",fill=NA),
      panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),panel.grid.minor=element_blank())
    print(gp)
    dev.off()
  }}

p.val.array <- matrix(nrow=3,ncol=3,dimnames=list(Metabolites,Comparison))
    
for (m in Metabolites){ 
    for (p in Comparison){
        sub.s <- subset(Corrected_Data, Metabolite==m & Group=="Susceptible")
        sub.r <- subset(Corrected_Data, Metabolite==m & Group=="Resistant")
        
        p.val.array[m,p] <- t.test(sub.s[[p]],sub.r[[p]])$p.val
    }}
    