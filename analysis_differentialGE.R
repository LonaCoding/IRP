library(limma)
library(reshape2)
library(dplyr)

BiocManager::install("sva")

library(sva)
library(ggplot2)
library(pheatmap)
library(biomaRt)

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

setwd("/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/")
rm(list=ls())
#ds=read.csv('~/Documents/VPA/Processed/VPA_data.csv')
ds=read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/Processed/VPA_data.csv')
#meta=read.csv('~/Documents/VPA/Info/GSM-Meta_combined_without_duplicated_GSMs_arrays_only.csv')
meta=read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/Info/GSM-Meta_combined_without_duplicated_GSMs_arrays_only.csv')

# preparing dataset for analysis
row.names(ds)=ds$ensembl_gene_id
ds2=ds[,-c(1:2)] #remove redundant ensembl_gene_id column
ds3=ds2[,colnames(ds2) %in% meta$GSM] #only select those listed in meta file
meta2=meta[,c(3,6,7,15:ncol(meta))]
row.names(meta2)=meta2$GSM
meta2=meta2[row.names(meta2) %in% colnames(ds3),]
meta3=meta2[order(row.names(meta2)),]
ds3=as.matrix(ds3[,order(colnames(ds3))])
identical(colnames(ds3),row.names(meta3))

# the reason for discrepancy in row numbers between ds3 and meta is that there was no raw data for GSE666

# some data visualisation
preCB_dens<-plotDensities(ds3)
preCB_mds<-plotMDS(ds3)
preCB_mds<-plotMDS(ds3, pch=16)
preCB_mds<-plotMDS(ds3, pch=16,gene.selection = "common")
preCB_bp<-boxplot(ds3)

meta3$batch=batch=paste(meta3$GPL,meta3$GSE,sep='_')
meta3$group2=(ifelse(meta3$group=='VPA',1,0))

ds4=normalizeBetweenArrays(ds3, method="quantile") #(Smyth and Speed, 2003): two-color microarray normalization -> limma
ds5=ComBat(ds4,batch = meta3$batch, mod = meta3$group2) #(Johnson, Li and Rabinovic, 2007) -> sva

# some data visualisation
posCB_dens<-plotDensities(ds5)
posCB_mds<-plotMDS(ds5)
posCB_mds<-plotMDS(ds5, pch=16)
posCB_mds<-plotMDS(ds5, pch=16,gene.selection = "common") #nornalise
posCB_bp<-boxplot(ds5)

# principal component analysis
pca <- prcomp(t(ds5))
as.data.frame(pca$x)->cmp # moving PC data into a separate object

ggplot(cmp,aes(x=cmp[,1],cmp[,2]))+geom_point(aes(colour=meta3$GSE))

# DE analysis
f <- factor(meta3$group)
design <- model.matrix(~0+f)
colnames(design) <- c("Control","VPA")

fit <- lmFit(ds5, design) #(Phipson et al., 2016) -> limma
contrast.matrix <- makeContrasts(VPA-Control,levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) #(Phipson et al., 2016) -> limma
summary(decideTests(fit2,adjust.method = "BH"))
topT=topTable(fit2, adjust="BH",number = nrow(ds3)) #(Phipson et al., 2016) -> limma
topT2=topT[topT$adj.P.Val<0.05,] # only get the diff expr ones that are sig

volcanoplot(fit2,highlight = 25)

pdf("volc1a.pdf",width = 5, height=5)
EnhancedVolcano(topT,lab=rownames(topT),x='logFC',y='adj.P.Val',shape = 9,pCutoff=0.05,ylab=bquote(~-Log[10] ~italic(FDR)), FCcutoff = 0,xlim=c(-2,2),axisLabSize=12,subtitle='',title='VPA vs Control',titleLabSize=12,legendPosition = 'none' , labSize=1, pointSize = 1, drawConnectors = FALSE )
dev.off()

pdf("volc1b.pdf",width = 5, height=5)
EnhancedVolcano(topT2,lab=rownames(topT2),x='logFC',y='adj.P.Val',pCutoff=0.05,ylab=bquote(~-Log[10] ~italic(FDR)), FCcutoff = 0,xlim=c(-2,2),axisLabSize=12,subtitle='',title='VPA vs Control',titleLabSize=12,legendPosition = 'none' , labSize=2 )
dev.off()

#plotlines()

write.csv(topT2,'~/Documents/GitRepo_Practiceals/IRP/VPA/a2/DE_analysis.csv')

write.csv(ds5,'~/Documents/GitRepo_Practiceals/IRP/VPA/a2/VPA_ComBated.csv')

#################################################################
#### gene sets analysis
#################################################################

read.csv('~/Documents/GitRepo_Practiceals/IRP/proc/Analysis/pathway annotations/Ensembl2Reactome_PE_All_Levels.txt',sep='\t',header = F)->pth
pth[pth$V8=='Homo sapiens',]->pth

read.csv('~/Documents/GitRepo_Practiceals/IRP/proc/Analysis/pathway annotations/GO_human_RNA_annotations.csv')->lnc_pth
read.csv('~/Documents/GitRepo_Practiceals/IRP/proc/Analysis/pathway annotations/RNAcentral_ensemblIDs.tsv',sep='\t',header = F)->rns_ids
data.frame('DB.Object.ID'=rns_ids$V1,
           'ensemblID'=do.call(rbind,strsplit(rns_ids$V6,'[.]'))[,1])->rns_ids2
unique(rns_ids2)->rns_ids2
do.call(rbind,strsplit(lnc_pth$DB.Object.ID,'_'))[,1]->lnc_pth$DB.Object.ID
merge(rns_ids2,lnc_pth,by='DB.Object.ID',all.y=T)->lnc_pth
lnc_pth[!is.na(lnc_pth$ensemblID),]->lnc_pth2

read.csv('~/Documents/GitRepo_Practiceals/IRP/proc/Analysis/pathway annotations/ChEBI2Reactome_PE_All_Levels.txt',sep='\t',header = F)->pthCh
pthCh[pthCh$V8=='Homo sapiens',]->pthCh

split(lnc_pth2,lnc_pth2$GO.ID)->lnc_pthL
lapply(lnc_pthL, function(x){x$ensemblID})->lnc_pthL2

split(pth,pth$V6)->pthL
lapply(pthL,function(x){x$V1})->pthL2

c(pthL2,lnc_pthL2)->pthL3

as.data.frame(ds5)->efitt
do.call(rbind,strsplit(row.names(efitt),'[.]'))[,1]->efitt$genes
lapply(pthL3,function(x){efitt$genes %in% x})->indices

# camera with not pre-ranked data (not the eBayes output)
camera(ds5,design,contrast= contrast.matrix,index=indices)->vpaAll
vpaAll2<-vpaAll
vpaAll2$path<-rownames(vpaAll2)
rownames(vpaAll2)<-NULL

ggplot(data=vpaAll2,aes(x=path,y=Pvalue,fill=path))+geom_bar() #vpaAll2 need to be dataframe
       
vpaAll[vpaAll$FDR<0.05 & complete.cases(vpaAll),]->vpaSig
write.csv(vpaSig,'~/Documents/GitRepo_Practiceals/IRP/VPA/a2/VPA_gene_set_enrichment.csv')
read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/a2/VPA_gene_set_enrichment.csv')->vpaSig
vpaSig2<-vpaSig
row.names(vpaSig2)->vpaSig2$Pathway

pthCh[pthCh$V6 %in% row.names(vpaSig2),]->vpa_Ch
pthL3[names(pthL3) %in% row.names(vpaSig2)]->vpaPths

lapply(vpaPths,function(x){x[x %in% row.names(topT)]->g
  if(length(g)==0){data.frame('gene_id'='gene','pth'=1)->h} else{data.frame('gene_id'=g,'pth'=1)->h}
  unique(h)
})->vpaGenes
vpaGenes2<-data.frame('gene_id'='gene','pth'=1)
for(n in 1:length(vpaGenes)){
  merge(vpaGenes2,vpaGenes[[n]],by='gene_id',all=T)->vpaGenes2
}
vpaGenes2[,-2]->vpaGenes2
colnames(vpaGenes2)[2:ncol(vpaGenes2)]<-names(vpaGenes)
unique(vpaGenes2)->vpaGenes2
topT$gene_id=row.names(topT)
merge(topT,vpaGenes2,by='gene_id',all.y=T)->vpaGenes3

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
geneNames=getBM(attributes = c('ensembl_gene_id','external_gene_name','description'),
                filters = 'ensembl_gene_id',
                values = vpaGenes3$gene_id, 
                mart = ensembl)
colnames(geneNames)[1]='gene_id'
vpaGenes3=merge(geneNames,vpaGenes3,by='gene_id',all.y=T)

write.csv(vpaGenes3,'~/Documents/GitRepo_Practiceals/IRP/VPA/a2/vpa_significant+pathway.csv',row.names = F)

row.names(vpaGenes2)=vpaGenes2$gene_id
vpaGenesH=as.matrix(vpaGenes2[,-1])
vpaGenesH2=apply(vpaGenesH,2,function(x){ifelse(is.na(x),0,x)})


pdf("vpaPathHMorig.pdf",width = 5, height=5)
origHm<-pheatmap(vpaGenesH2) #original. row label need fixing
dev.off()


vpaGenesH2_DF <- as.data.frame(vpaGenesH2)
vpaGenesH2_DF$gene_id<-rownames(vpaGenesH2_DF)
ensembl2 <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
geneNames_vpaGenesH2=getBM(attributes = c('ensembl_gene_id','external_gene_name','description'),
                filters = 'ensembl_gene_id',
                values = vpaGenesH2_DF$gene_id, 
                mart = ensembl2)
colnames(geneNames_vpaGenesH2)[1]='gene_id'
vpaGenesH2m<-merge(vpaGenesH2_DF,geneNames_vpaGenesH2,by='gene_id',all.x=T,all.y=F)

rownames(vpaGenesH2m) = make.names(vpaGenesH2m$external_gene_name, unique=TRUE)
vpaGenesH2m_gl<-vpaGenesH2m$external_gene_name
vpaGenesH2m_glrn<-rownames(vpaGenesH2m)
end<-length(vpaGenesH2m)-2
vpaGenesH2mc<-vpaGenesH2m[1:end]
vpaGenesH2mc2<-vpaGenesH2mc[-c(1)]


modHm<-pheatmap(vpaGenesH2mc2,fontsize_row=3,fontsize_col =3.5,filename="vpaPathwayHeatmapGeneName.png")#,scale = "row",filename="vpaPathwayHeatmap.png")
pdf("vpaPathHM.pdf",width = 5, height=5)
pheatmap(vpaGenesH2mc2,fontsize_row=3,fontsize_col =3.5)#,filename="vpaPathwayHeatmapGeneName.png")#,scale = "row",filename="vpaPathwayHeatmap.png")
dev.off()


