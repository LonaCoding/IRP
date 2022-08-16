################################
# extracting NimbleGen data
################################

setwd('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/rawDataNI/nimblegen/GSE65297/') # NimbleGen data after untarring

## The section below, upto line 34, is entirely written by Dr. Marcin Wozniak, detailing how to manually create an annotation package for illumina: ####
# there is no annotation package thaht can be downloaded automatically from Bioconductor
# need to build and annotation package with 'pdInfoBuilder'
# here is info on how to do it: 
# https://www.bioconductor.org/packages/release/bioc/vignettes/pdInfoBuilder/inst/doc/BuildingPDInfoPkgs.pdf

#BiocManager::install("pdInfoBuilder",force = TRUE)
#
#library(pdInfoBuilder)
#ndf=list.files('.',pattern = 'ndf$') # I think the ndf file needs to be unzipped
#xys=list.files('.',pattern = 'xys$') # the same here, make sure that there is only 
#                                     # 1 item in the xys vector, otherwise you're gonna get an error
#
#seed <- new("NgsExpressionPDInfoPkgSeed",
#            ndfFile = ndf, 
#            xysFile = xys,
#            author = "Marcin Wozniak",
#            email = "mw299@le.ac.uk",
#            biocViews = "AnnotationData",
#            organism = "Human", 
#            species = "Homo Sapiens")
#
#makePdInfoPackage(seed, destDir = ".")

# Now need to install the package 
#install.packages('./pd.gpl16025.100718.hg18.opt.expr/',repos = NULL) #make sure dir is correct (to where package is)

######

install.packages('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/pd.gpl16025.100718.hg18.opt.expr/',repos = NULL) #make sure dir is correct (to where package is)



install.packages("BiocManager")

BiocManager::install("oligo")
BiocManager::install("limma")
BiocManager::install("beadarray")
BiocManager::install("lumi")
BiocManager::install("GEOquery")

library(oligo)
library(limma)
library(beadarray)
library(lumi) 
library(R.utils)
library(GEOquery) 

# all the above needs to be done only once
###########################################################

# NimbleGen data comes in xys files (no need to gunzip them)
# they can be read with oligo package

fls=list.files('.',pattern = 'xys')
z=read.xysfiles(fls,pkgname='pd.gpl16025.100718.hg18.opt.expr')
y=rma(z,normalize=F)

#write.exprs(y,'~/Documents/VPA/Processed/GSE65297[100718_HG18_opt_expr].csv')
write.exprs(y,'/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/Processed/GSE65297[100718_HG18_opt_expr].csv')

################################
# extracting Illumina data
################################

setwd('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/rawDataNI/illumina/GSE60154/') 
# The first 4 lines from the GSE60154_non-normalized.txt file need to be removed leving only the expression data
x<-read.ilmn('GSE60154_non-normalized_mod.txt',expr = '14P0100390',annotation = 'ID_REF')
nec(x)->y # in-chip normalisation
log2(y$E)->E # log2 transformation of expression data
data.frame(y$genes,E)->z # adding Illumina probe IDs
#write.csv('~/Documents/VPA/Processed/GSE60154_[HT-12 V4.0].csvz')
#write.csv('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/Processed/GSE60154_[HT-12 V4.0].csvz')
write.table(z,'/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/Processed/GSE60154_[HT-12 V4.0].csv',se='\t')

########

setwd('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/rawDataNI/illumina/GSE59803/') 
# The first 4 lines from the GSE60154_non-normalized.txt file need to be removed leving only the expression data
#x<-read.ilmn('GSE59803_non-normalized_mod.txt',expr = '14P0100390',annotation = 'ID_REF')
x<-read.ilmn('GSE59803_non-normalized_mod.txt',expr = 'SAMPLE',annotation = 'ID_REF')
nec(x)->y # in-chip normalisation
log2(y$E)->E # log2 transformation of expression data
data.frame(y$genes,E)->z # adding Illumina probe IDs
#write.csv('~/Documents/VPA/Processed/GSE60154_[HT-12 V4.0].csvz')
#write.csv('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/Processed/GSE59803_[HT-12 V4.0].csvz')
write.table(z,'/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/Processed/GSE59803_[HT-12 V4.0].csv',se='\t')
