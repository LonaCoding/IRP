library(dplyr)
library(ggplot2)
library(pheatmap)

setwd("/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/")
#rm(list=ls())

#vPath=read.csv('~/Documents/VPA/analysis/vpa_significant+pathway.csv')
#oPath=read.csv('~/Documents/VPA/analysis/DEgenes_Ob35.OverOb.csv')
#nPath=read.csv('~/Documents/VPA/analysis/DEgenes_OverOb.UnderNorm.csv')

vPath=read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/analysis/vpa_significant+pathway.csv')
oPath=read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/analysis/DEgenes_Ob35.OverOb.csv')
nPath=read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/analysis/DEgenes_OverOb.UnderNorm.csv')

#vErch=read.csv('~/Documents/VPA/analysis/VPA_gene_set_enrichment.csv')
#oErch=read.csv('~/Documents/VPA/analysis/pth_enrichment_OverOb.UnderNorm.csv')
#nErch=read.csv('~/Documents/VPA/analysis/pth_enrichment_Ob35.OverOb.csv')

vErch=read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/analysis/VPA_gene_set_enrichment.csv')
oErch=read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/analysis/pth_enrichment_OverOb.UnderNorm.csv')
nErch=read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/analysis/pth_enrichment_Ob35.OverOb.csv')


## pathway

makePathHeatmap <- function (a,b,al,bl){
  t=c(al,bl)
  
  #comparing VPA with class2-vs-class1overweight, then only getting the GENES that are present in both datasets and merging column-wise
  ab<-merge(a,b, by.x = "gene_id", by.y = "ensembl_gene_id")
  
  #removing duplicate columns. 
  abCln<-subset(ab, select=-c(external_gene_name.y,description.y))
  
  #rearrange such that logFCs are grouped together and so on
  abGrp<-abCln%>%relocate(logFC.y, .after=logFC.x)%>%relocate(AveExpr.y, .after=AveExpr.x)%>%
    relocate(t.y, .after=t.x)%>%relocate(P.Value.y, .after=P.Value.x)%>%
    relocate(adj.P.Val.y, .after=adj.P.Val.x)%>%relocate(B.y, .after=B.x)%>%relocate(hgnc_id, .after=B.y)
  
  #get only logFC. P value not needed because asssumed everything there is p<0.05 or significant
  abReg<-abGrp[,c(1:9,14,15)]
  abFC<-abGrp[,c(1:5)]
  abFCFinal=abFC
  
  #do comparison between fold changes logFCs. Find consensus (both are upreg or downreg at same time). put in new column
  #mutate(abFCFinal, consensus=compareFC(logFC.x,logFC.y))
  #apply(abFCFinal,1, compareFC, abFCFinal$logFC.x,abFCFinal$logFC.y)
  
  
  abFCb<-abFC[,c(2,4,5)]
  row.names(abFCb)=abFC$external_gene_name.x
  abFCmat<-as.matrix(abFCb[,c(2,3)])
  print(t)
  h<-pheatmap(abFCmat, display_numbers = T, labels_col=c(al, bl), angle_col=0)
  #font color based on group?
  return (abFCFinal)
}

makeErchHeatmap <- function (a,b,al,bl,outType){
  t=c(al,bl)
  
  #comparing VPA with class2-vs-class1overweight, then only getting the GENES that are present in both datasets and merging column-wise
  ab<-merge(a,b, by.x = "X", by.y = "X")
  
  #rearrange such that logFCs are grouped together and so on
  abGrp<-ab%>%relocate(Direction.x, .after=X)%>%relocate(Direction.y, .after=Direction.x)%>%relocate(NGenes.x, .after=Direction.y)%>%relocate(NGenes.y, .after=NGenes.x)%>%relocate(FDR.x, .after=NGenes.y)%>%relocate(FDR.y, .after=FDR.x)%>%
    relocate(PValue.x, .after=FDR.y)%>%relocate(PValue.y, .after=PValue.x)
  
  #get only logFC. P value not needed because asssumed everything there is p<0.05 or significant
  abFDR<-abGrp[,c(1:7)]
  abShort<-abGrp[,c(1:5)]
  
  #do comparison between fold changes logFCs. Find consensus (both are upreg or downreg at same time). put in new column
  #mutate(abFCFinal, consensus=compareFC(logFC.x,logFC.y))
  #apply(abFCFinal,1, compareFC, abFCFinal$logFC.x,abFCFinal$logFC.y)
  
  
  #abFCb<-abFC[,c(2,4,5)]
  abFC<-abGrp
  abFCb<-abGrp
  
  row.names(abFCb)=abFC$X
  abFCbBin<-abFCb
  colnames(abFCbBin)[2]<-"Direction_Name.x"
  colnames(abFCbBin)[3]<-"Direction_Name.y"
  abFCbBin$Direction.x[abFCbBin$Direction_Name.x=='Up']<-1
  abFCbBin$Direction.x[abFCbBin$Direction_Name.x=='Down']<--1
  abFCbBin$Direction.y[abFCbBin$Direction_Name.y=='Up']<-1
  abFCbBin$Direction.y[abFCbBin$Direction_Name.y=='Down']<--1
  abFCbBin2<-abFCbBin%>%relocate(Direction.x, .after=Direction_Name.y)%>%relocate(Direction.y, .after=Direction.x)
  abFCmatNG<-as.matrix(abFCbBin2[,c(6,7)]) #Number of Genes
  abFCmatDir<-as.matrix(abFCbBin2[,c(4,5)]) #Direction
  
  if (outType==0){
    h<-pheatmap(abFCmatNG, display_numbers = T, labels_col=c(al, bl),angle_col=0)
  } else {
    h2<-pheatmap(abFCmatDir, display_numbers = T, labels_col=c(al, bl),angle_col=0) 
  }
   
  return (abFC)
}

voFCres<-makePathHeatmap(vPath,oPath,"VPA","BMI >= 35 vs 25<=BMI<35")
vnFCres=makePathHeatmap(vPath,nPath,"VPA","25<=BMI<35 vs BMI<25") #REPRINT!!!

vnFCresNG<-makeErchHeatmap(vErch,nErch,"VPA","25<=BMI<35 vs BMI<25",0) #Number of Genes
vnFCresD<-makeErchHeatmap(vErch,nErch,"VPA","25<=BMI<35 vs BMI<25",1) #Direction
voFCresNG<-makeErchHeatmap(vErch,oErch,"VPA","BMI >= 35 vs 25<=BMI<35",0) #Number of Genes
voFCresD<-makeErchHeatmap(vErch,oErch,"VPA","BMI >= 35 vs 25<=BMI<35",1) #Direction

#Direction:
#1 up
#-1 down

# positive logfc means upreg in o where Ob35 (BMI>25) (the first) while vs overweight, then in n is overweight while vs underweight 