library(tcltk)
library(ggplot2)

file <- tk_choose.files(default="/home/sdenecke/Dropbox/PhD_Projects/DGRP_Paper/Results/Cyp6g1KO/",multi=F)
Data <- read.csv(file)
setwd(dirname(file))
dir.create("./GraphDir")
setwd("./GraphDir")

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


Metabolites=unique(Corrected_Data$Metabolite)
Comparison=c("Media","Larvae")
Corrected_Data <- subset(Corrected_Data,set=="RAL517")

for (m in Metabolites){
    for (p in Comparison){
        sub=subset(Corrected_Data,Metabolite==m)
        lim <- ceiling(max(sub[["Larvae"]],sub[["Media"]])*1.2)
        pdf(file=paste(m,p,"RAL517","pdf",sep="."))
        
        gp=ggplot(sub,aes_string(x='Genotype',y=p,fill="Genotype")) 
        gp=gp+geom_dotplot(binaxis='y', stackdir='center') 
        gp=gp+ylab(paste("Concentration of ",m," [ppb]",sep=""))
        gp=gp+stat_summary(geom="errorbar",colour='black',size=1,width=.4,fun.data="mean_se")
        gp=gp+stat_summary(fun.y="mean",fun.ymax="mean",fun.ymin="mean", geom="crossbar", width=.4,colour='black',size=.8)
        gp=gp+ylim(c(0,lim))
        gp=gp+ggtitle(paste(m,p))
        gp=gp+scale_fill_manual(values=c("grey95","grey5"))
        gp=gp+theme(text=element_text(size=18,family="",face="bold"),
                    axis.text.x=element_text(size=16,colour="black"),
                    axis.text.y=element_text(size=14,family="",face="bold"),
                    panel.background=element_rect(fill=NA),
                    panel.border=element_rect(colour="black",fill=NA),
                    panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),legend.position="none")
        print(gp)
        dev.off()
    }}

p.val.array <- matrix(nrow=3,ncol=3,dimnames=list(Metabolites,Comparison))

for (m in Metabolites){ 
    for (p in Comparison){
        
        sub.s <- subset(Corrected_Data, Metabolite==m & Genotype=="Wxac")
        sub.r <- subset(Corrected_Data, Metabolite==m & Genotype=="Wxac[Cyp6g1KO]")
        
        p.val.array[m,p] <- t.test(sub.s[[p]],sub.r[[p]])$p.val
    }}

