#Import Libraries
library(plyr)
library(lme4)
library(RColorBrewer)


# Import your data from .csv file
File_Name=file.choose()
mainDir=dirname(File_Name)
setwd(mainDir)
Data=read.csv(File_Name)

###Let's get a workable Data frame from your raw Data File
{         
SPLIT_Frame=matrix(unlist(strsplit(as.character(Data$Image.Name),"\\).")),ncol=8,byrow=T)[,c(1:6,8)]
SPLIT_Frame=gsub("\\(","",SPLIT_Frame)
SPLIT_Frame=gsub(".tif","",SPLIT_Frame)
SPLIT_Frame=gsub("opRight","TopRight",SPLIT_Frame)
SPLIT_Frame=gsub("opLeft","TopLeft",SPLIT_Frame)
SPLIT_Frame=gsub("owRight","LowRight",SPLIT_Frame)
SPLIT_Frame=gsub("owLeft","LowLeft",SPLIT_Frame)
WI_Frame=cbind(SPLIT_Frame,Data$Wiggle.Index)
colnames(WI_Frame)=c("Genotype","Dose","Time","Plate","Insecticide","Date","Position","Wiggle_Index")
WI_Frame[,'Time']=as.numeric(gsub('min',"",WI_Frame[,'Time']))		  
WI_Frame=data.frame(WI_Frame) 
}
    
    
    ## Ok now lets run this horribly inefficient way to get corrected values for each genotype
    {
    
    Corrected_Data=data.frame()
    Dates=as.character(unique(WI_Frame$Date))
    Doses=as.character(unique(WI_Frame$Dose))
    Genotypes=as.character(unique(WI_Frame$Genotype))
    Plates=as.character(unique(WI_Frame$Plate))
    Positions=as.character(unique(WI_Frame$Position))
    Insecticides=as.character(unique(WI_Frame$Insecticide))
    for (a in Dates){
      for (b in Doses){
        for(c in Genotypes){
          for (d in Plates){
            for (e in Positions){
              for (f in Insecticides){
                sub=subset(WI_Frame,Date==a & Dose==b & Genotype==c & Plate==d & Position==e & Insecticide==f)
                if(length(sub$Wiggle_Index)>0){
                  subWI=as.numeric(as.character(sub$Wiggle_Index))
                  EV=c()
                  for (i in 1:length(subWI)){
                    x=subWI[i]/subWI[1]
                    EV[i]=x
                  }
                  sub$Modified=EV
                  Corrected_Data=rbind(Corrected_Data,sub)
                }
              }
            }
          }
        }
      }
    }
    }

## Lets save this data frame so you can call it later if need be. The previous analysis can take a while
write.csv(Corrected_Data,file="Completed_Data_Frame.csv")
Corrected_Data=read.csv("Completed_Data_Frame.csv")


# You might want to check your data frame to see that nothing strange is happening (e.g duplicated values)
table(Corrected_Data$Genotype)[which(table(Corrected_Data$Genotype)!=32)]


## Ok now lets see if transforming your data helps its normalization at all
Corrected_Data=subset(Corrected_Data,Time==60)
Dose=unique(Corrected_Data$Dose)
Transformation=c(raw="",log="log",sqrt="sqrt")
norm.test=list()
Vg=list()
 BLM.l=list()
 BLUP.l=list()
for (a in Dose){
  norm.test[[a]]=list()
  Vg[[a]]=list()
  Corrected_Data.Dose=subset(Corrected_Data,Dose==a)
    for (c in Transformation) {
       if (c=="") Data.Trans=Corrected_Data.Dose[,'Modified']
      if (c=="log") Data.Trans=log(Corrected_Data.Dose[,'Modified'])
      if (c=="sqrt") Data.Trans=sqrt(Corrected_Data.Dose[,'Modified'])
      norm.test[[a]][[names(Transformation)[Transformation==c]]]=as.numeric(shapiro.test(Data.Trans)[2])
      test=lmer(Data.Trans~1+(1|Corrected_Data.Dose$Genotype))
      Vg[[a]][[names(Transformation)[Transformation==c]]]=c(Vg=NA,Ve=NA)
      Vg[[a]][[names(Transformation)[Transformation==c]]]["Vg"]=VarCorr(test)$'Corrected_Data.Dose$Genotype'[1]
      Vg[[a]][[names(Transformation)[Transformation==c]]]["Ve"]=Vg[[a]][[names(Transformation)[Transformation==c]]]["Vg"]+attr(VarCorr(test),"sc")^2
    } 
  }

# Lets Calculate heritability
H=list()
  for (a in Dose){
    H[[a]]=list()
      for (c in Transformation) {
        r=names(Transformation)[Transformation==c]
        gen=Vg[[a]][[r]]
        H[[a]][[r]]['H']=gen[['Ve']]/(gen[['Ve']]+gen[['Vg']])
}
    }
write.csv(unlist(H),file="heritability.csv")

# Let's calculate the best linear unbiased predictor or BLUP  
for (a in Dose){
  Corrected_Data.Dose=subset(Corrected_Data,Dose==a)
  BLM=BLM=lmer(Modified~1+(1|Genotype),data=Corrected_Data.Dose)
  BLM.l[[a]]=BLM
  shapiro.test(ranef(BLM.l[[a]])$Genotype[,1])
  Genotype=as.character(sort(unique(Corrected_Data.Dose$Genotype)))
  Genotype=gsub("RAL","",Genotype)
  Genotype=as.numeric(Genotype)
  BLUP=(summary(BLM.l[[a]])$coefficients[1,1]+ranef(BLM.l[[a]])$Genotype[,1])^2
  Polished=na.omit(data.frame(cbind(Genotype,BLUP)))
  BLUP.l[[a]]=Polished
  #write.csv(Polished,row.names=F,file=paste('BLUP',a,'.csv'))
}  

l=list(Doses)
for (a in Dose){
Corrected_Data.Dose=subset(Corrected_Data,Dose==a)
Summary.total=ddply(Corrected_Data.Dose, c("Genotype","Time"), summarise,
            N    = length(Modified),
            mean = mean(Modified),
            sd   = sd(Modified),
            se   = sd / sqrt(N),
            ci   = qnorm(0.975)*se)
l[[a]]=Summary.total$mean
Summary.total$BLUP=BLUP.l[[a]]$BLUP
Summary.total=Summary.total[order(Summary.total$mean),]
Color=colorRampPalette(c("steelblue","orangered"))
Summary.total$rc <-  Color(length(Summary.total$Genotype))
#Summary.total$rc=rev(grey((Summary.total$mean/max(Summary.total$mean))))

par(mar=c(5,4,2,1),mgp=c(2.2,.5,0),mfrow=c(1,1))
pdf(file=paste("Range of DGRP Response at",a,"Imidacloprid.pdf"),width=10)
bp=barplot(Summary.total$mean,col=Summary.total$rc,ylim=c(0,.05+max(Summary.total$mean+Summary.total$ci)),border=TRUE,axes=F,family="serif",ylab=list("RMR Value",cex=1.5,font=2),main=a,cex.main=2)
segments(bp,Summary.total$mean-Summary.total$ci, bp,Summary.total$mean+Summary.total$ci,col="black",lwd=2)
segments(x0=bp-.3, y0=Summary.total$mean+Summary.total$ci,x1=bp+.3,lwd=2,col="black")
axis(side=2,at=c(0,.2,.4,.6,.8,1.0),family="serif",font=2,las=1,tick=T)
mtext("DGRP Genotype",side=1,cex=1.5,font=2,family="serif",line=1)
#grid(nx=NA,ny=NULL)
par(new=T)
bp=barplot(Summary.total$mean,col=Summary.total$rc,ylim=c(0,.05+max(Summary.total$mean+Summary.total$ci)),axes=F,family="serif")
dev.off()

write.csv(data.frame(Summary.total$Genotype,Summary.total$BLUP),file=paste("BLUP GWAS Input Values for",a,".csv"))
write.csv(data.frame(Summary.total$Genotype,Summary.total$mean),file=paste("Mean GWAS Input Values for",a,".csv"))



} 

cor.spearman <- cor(l$'25ppm',l$'100ppm',method="spearman")
Corr_Spearman <- cor.test(l$'25ppm',l$'100ppm',method="spearman")
Corr_Kendall <- cor.test(l$'25ppm',l$'100ppm',method="kendall")
Corr_Pearson <- cor.test(l$'25ppm',l$'100ppm',method="pearson")


