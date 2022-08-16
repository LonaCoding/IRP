#weighted correlation network analysis 




#library(R.utils)


## WGCNA install preparation

#install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "BiocManager")
#BiocManager::install(c("GO.db", "preprocessCore", "impute")); 

##old ver in case bioconductor fails
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GO.db", "preprocessCore", "impute")) 

## if going to use annotation stuff
#orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg");
#orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6));
#packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep="");
#
#BiocManager::install(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList")) 


# WGCNA install
install.packages("BiocManager");
install.packages("fastDummies");
BiocManager::install("WGCNA");

library(fastDummies); #needed to one-hot encode the group (vpa/control)
library(ggplot2)
library(reshape2)

# manual download-install (just in case)
#install.packages("path/to/file", repos = NULL, lib = "path/to/library")

#run on every new seesion
library(WGCNA);
setwd('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/')
options(stringsAsFactors = FALSE); #see tutorial.

## data preparation stage



# gene expression data set, with duplicate and batch effect removed.
#data must be: one gene per row, one GSE per column
inGEdata = read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/Processed/VPA_ComBated.csv');

#basic info of input dataset
dim(inGEdata);
names(inGEdata);

# store data
GEdf = as.data.frame(t(inGEdata[, -c(1)])); #transpose (switch columns and rows) while removing gene IDs
names(GEdf) = inGEdata$X; #get list of genes IDs and assign them as column names
rownames(GEdf) = names(inGEdata)[-c(1)];




# check if data is ready
gsg = goodSamplesGenes(GEdf, verbose = 3);
gsg$allOK #if true, says TRUE. if not, removes problematic data.

#only use this if above code says FALSE
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(GEdf)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(GEdf)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  GEdf = GEdf[gsg$goodSamples, gsg$goodGenes]
}

# clustering to plot dendogram to find outliers
sampleTree = hclust(dist(GEdf), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# remove outlier (only run if any found) # justify why doing this.
# Plot a line to show the cut
abline(h = 60, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
finalGEdf = GEdf[keepSamples, ]
nGenes = ncol(finalGEdf)
nSamples = nrow(finalGEdf)

#the outlier removal code isn't run because there are no outliers. Even those above 60 in the plot may be considered to be not outliers.


#  clinical data=meta csv file
# clinical data file must contain: VPA dosage and VPA/control for each GSE (1 GSE per row)
inCTdata= read.csv('~/Documents/GitRepo_Practiceals/IRP/VPA/Info/GSM-Meta_combined_without_duplicated_GSMs_arrays_only.csv');

#basic info of clinical train input dataset
dim(inCTdata)
names(inCTdata)

# remove columns that hold information we do not need.
CTonly = inCTdata[, c(7,15) ]; #only take GSM and VPA?control. VPA concentration can be done later once it has been standardized into numeric.
#CTonlyFinal<-cbind(to.dummy(CTonly$group, "Group")) #alt method, need to append column using cbind and delete original group column. need to find package

CTonlyDum<-fastDummies::dummy_cols(CTonly,select_columns="group") #one-hot encode group column to turn it into numeric
CTonlyFinal=CTonlyDum[,-c(2)] #data only contains GSM and binary group columns

dim(CTonly)
names(CTonly)

dim(CTonlyFinal)
names(CTonlyFinal)

# Form a data frame analogous to expression data that will hold the clinical traits.
finalGEdf=GEdf; #only do this if no outlier removal was done
sample = rownames(finalGEdf);
traitRows = match(sample, CTonlyFinal$GSM); #GSM names
datTraits = CTonlyFinal[traitRows, -1];
rownames(datTraits) = CTonlyFinal[traitRows, 1]; #final

collectGarbage();


#see prep dendogram
# Re-cluster samples
sampleTree2 = hclust(dist(finalGEdf), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

save(finalGEdf, datTraits, file = "prepWGCNA.RData") #backuping work so far

## actual network making

#enableWGCNAThreads() # skip if running RStudio 
disableWGCNAThreads()

prep = load(file ="prepWGCNA.RData") #load saved work

DE=read.csv('~/Documents/GitRepo_Practiceals/IRP/DE_analysis.csv')

finalGEdf2=finalGEdf[,colnames(finalGEdf) %in% DE$X]


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 
#explain why use this parameter value (default: from = 12, to=20, by=2)
#these are possible powers that can be tried to see which one is best. 

# Call the network topology analysis function (this bit kept making Rstudio crash)
sft = pickSoftThreshold(finalGEdf2, powerVector = powers, verbose = 5) #explain why use this parameter value (default: verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# make network
chosenPower=sft$powerEstimate #set the power according to results above (refer to network topology graph, find the best fit, in which its lowest value where the graph curve flattens out)
net = blockwiseModules(finalGEdf, power = chosenPower, #explain parameter value choices
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "vpaTOM", 
                       verbose = 3)

#default params
#power=6 # very dependent on network topology graph. avoid using this value.
#TOMType = "unsigned", 
#minModuleSize = 30,
#reassignThreshold = 0, 
#mergeCutHeight = 0.25,
#numericLabels = TRUE, 
#pamRespectsDendro = FALSE,
#saveTOMs = TRUE,
#saveTOMFileBase = "vpaTOM", #adapted
#verbose = 3)


## display output
## not needed
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

## merging networks
adjacency = adjacency(finalGEdf2, power = chosenPower);
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average"); #can be loaded from networkOutWGCNA.Rdata

prep2 = load(file ="networkOutWGCNA.RData") #load saved work
prep3 = load(file ="prepWGCNA.RData") #load saved work


plot(geneTree)

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE);
table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(finalGEdf2, colors = dynamicColors)
MEs = MEList$eigengenes #can be loaded from networkOutWGCNA.Rdata



# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");#can be loaded from networkOutWGCNA.Rdata

plot(METree)

MEDissThres = 0.2

abline(h=MEDissThres, col = "red")

merge = mergeCloseModules(finalGEdf2, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs; #final

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)



#######

identical(row.names(datTraits), row.names(MEs)) 


datTraits$ID=row.names(datTraits)

MEs$ID=row.names(MEs)

MEsdatTraits<-merge(datTraits,MEs,by='ID',all=T)

MEsdatTraits$VPA=ifelse(MEsdatTraits$group_VPA==1,'VPA','CTRL')


MEsdatTraits2=MEsdatTraits[,-c(2,3)]


MEsdatTraits3=melt(MEsdatTraits2,by=c('ID','VPA'))

ggplot(MEsdatTraits3,aes(x=VPA,y=value))+geom_boxplot()+facet_wrap(~variable,scales='free')

#boxplot(MEsdatTraits2)

## for 
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


cbind('Network'=mergedColors,'X'=colnames(datExpr))->vGenes
ranks2=merge(ranks,vGenes,by='X',all=T)
write.csv(ranks2,'Analysis/DE_analysis+Nets.csv',row.names = F)

# plotting eigengene values
mergedMEs$GSM=row.names(mergedMEs)
MEs=merge(meta,mergedMEs,by='GSM')

melt(MEs,id=c("GSM","group"))->MEsplotM

ggplot(MEsplotM,aes(x=group,y=value))+
  geom_jitter(aes(colour=group),width=0.2)+
  geom_boxplot(aes(fill=group))+
  facet_wrap(~variable)


#shapiro.test(log(MEsdatTraits$MEred)) #if p<0.05, data not normally distributed

#use wilcoxon test on eigengenevalues if p<0.05=data not normally distributed. else t.test (same syntax as wilcox)





#do this fore each networks/MEs
#wRed<-wilcox.test(MEsdatTraits[MEsdatTraits$VPA=='VPA',]$MEred,MEsdatTraits[MEsdatTraits$VPA=='CTRL',]$MEred)

#old guideline:
#col 1= nuetwork name / MEs
#col 2= is it dist normally or not
#col 3= mean / median (depend on distr)
#col 4= type of test (tetst/wilcox)
#col 5= pvalue
#col 6= significant?

#get list of networks/MEs
listMEs=colnames(MEsdatTraits)[5:ncol(MEsdatTraits)-1]

# making the dataframe containing the significant result
significantResList<-data.frame(Network_Name=vector(mode="character",length = length(listMEs)), #1
                          Shapiro_Pvalue=vector(mode="double",length = length(listMEs)), #2
                          NormalDistributed=vector(mode="logical",length = length(listMEs)), #3: is it dist normally or not
                          Mean_Median_VPA=vector(mode="double",length = length(listMEs)), #4: mean if dist normally, else median
                          Mean_Median_Control=vector(mode="double",length = length(listMEs)), #5: mean if dist normally, else median
                          SD_IQR_VPA=vector(mode="double",length = length(listMEs)), #6: SD if dist normally, else IQR
                          SD_IQR_Control=vector(mode="double",length = length(listMEs)), #7: SD if dist normally, else IQR
                          Test_Type=vector(mode="logical",length = length(listMEs)), #8: use Wilcoxon or T-test
                          Test_Pvalue=vector(mode="double",length = length(listMEs)), #9: test result of T-test if dist normally, else Wilcoxon
                          SignificantResult=vector(mode="logical",length = length(listMEs))) #10: is it significant or not




i=1
for (n in listMEs){ #for each network, whereby n is the network name
  print(n)
  significantResList[i,1]=n #Network_Name
  
  # check wether they are normally distributed or not, using shapiro test
  shapiroRes<-shapiro.test(log(MEsdatTraits[[n]])) #if p<0.05, data not normally distributed
  shapiroPval<-shapiroRes$p.value
  significantResList[i,2]=shapiroPval #Shapiro_Pvalue

  
  if (shapiroPval<0.05){
    #NOT normally distributed
    significantResList[i,3]=FALSE #NormalDistributed
    
    #median and IQR if NOT standard distributed shapiro.test p<0.05
    significantResList[i,4]<-median(MEsdatTraits[MEsdatTraits$VPA=='VPA',][[n]]) #Mean_Median_VPA
    significantResList[i,5]<-median(MEsdatTraits[MEsdatTraits$VPA=='CTRL',][[n]]) #Mean_Median_Control
    significantResList[i,6]<-IQR(MEsdatTraits[MEsdatTraits$VPA=='VPA',][[n]]) #SD_IQR_VPA
    significantResList[i,7]<-IQR(MEsdatTraits[MEsdatTraits$VPA=='CTRL',][[n]]) #SD_IQR_Control
    
    significantResList[i,8]="Wilcoxon" #Test_Type
    testRes<-wilcox.test(MEsdatTraits[MEsdatTraits$VPA=='VPA',][[n]],MEsdatTraits[MEsdatTraits$VPA=='CTRL',][[n]])$p.value
    significantResList[i,9]<-testRes #Test_Pvalue
    if (testRes<0.05){
      #significant result
      significantResList[i,10]=TRUE #SignificantResult
    } else {
      #NOT significant result
      significantResList[i,10]=FALSE #SignificantResult
    }
    
    
  } else {
    # Normally distributed
    wilcoxResList[i,3]=TRUE #NormalDistributed
    
    #mean and sd if standard dist = shapiro.test p>0.05
    significantResList[i,4]<-mean(MEsdatTraits[MEsdatTraits$VPA=='VPA',][[n]]) #Mean_Median_VPA
    significantResList[i,5]<-mean(MEsdatTraits[MEsdatTraits$VPA=='CTRL',][[n]]) #Mean_Median_Control
    significantResList[i,6]<-sd(MEsdatTraits[MEsdatTraits$VPA=='VPA',][[n]]) #SD_IQR_VPA
    significantResList[i,7]<-sd(MEsdatTraits[MEsdatTraits$VPA=='CTRL',][[n]]) #SD_IQR_Control
    
    significantResList[i,8]="T" #Test_Type
    testRes<-t.test(MEsdatTraits[MEsdatTraits$VPA=='VPA',][[n]],MEsdatTraits[MEsdatTraits$VPA=='CTRL',][[n]])$p.value
    significantResList[i,9]<-testRes #Test_Pvalue
    if (testRes<0.05){ #p<0.05 means significant
      #significant result
      significantResList[i,10]=TRUE #SignificantResult
    } else {
      #NOT significant result
      significantResList[i,10]=FALSE #SignificantResult
    }
    
  }
  
  #MEsdatTraits[MEsdatTraits$VPA=='VPA',][n]

  i=i+1
}

significantResList2<-significantResList

prep4 = load(file ="sigResAnalysisWGCNA.RData") #load saved work

save(significantResList, 
     file = "sigResAnalysisWGCNA.RData")


#p<0.05 means significant
## for each MEs, from eigenvalues, for VPA
#mean and sd if standard dist = shapiro.test p>0.05

#median and IQR if NOT standard distributed shapiro.test p<0.05


##using wilcox test p value result, for those that are siginificant diff (p<0.05), 
##in cytoscape. install stringAPp. copy paste list of ENSMBL genes (if fail, go to biomart and convert)

## to do:
#make a table
#col 1= nuetwork name / MEs
#col 2= is it dist normally or not
#col 3= mean / median (depend on distr)
#col 4= type of test (tetst/wilcox)
#col 5= pvalue
#col 6= significant?



#dynamicColors
DE2<-DE
geneNet=data.frame('X'=colnames(finalGEdf2), 'Network'=dynamicColors); merge(geneNet,DE,by='ID',all=T)
DE2$NET=dynamicColors

#cytoscape network. 
# identify genes that is in networks. use DE2
# for each network, get list of genes from that and put through cytoscape




## save data

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "networkOutWGCNA.RData")

save(DE, DE2, finalGEdf2, finalGEdf, TOM, dissTOM, geneTree, geneNet, mergedMEs, MEs, MEDiss,MEDissThres,datTraits,METree, MEsdatTraits, MEsdatTraits3,
     file = "networkOutWGCNADEdatTraits.RData")

save(DE, DE2, finalGEdf2, finalGEdf, geneTree, geneNet,
     file = "networkOutWGCNA_DEgeneTree.RData")

save(TOM, dissTOM,
     file = "networkOutWGCNA_TOM.RData")

save(MEs, MEDiss,MEDissThres,datTraits,METree, MEsdatTraits, MEsdatTraits3,
     file = "networkOutWGCNAMEdatTraits.RData")

prep2 = load(file ="networkOutWGCNA.RData") #load saved work
#MEs = eigen values


brownNet=read.csv('../networkSig/brown.csv')
brownNet=read.csv('../networkSig/brown.csv')
brownNet=read.csv('../networkSig/brown.csv')
brownNet=read.csv('../networkSig/brown.csv')


