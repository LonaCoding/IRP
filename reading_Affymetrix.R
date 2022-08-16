library(limma)
library(R.utils)
library(affy)
library(oligo)

rm(list=ls())
#setwd('~/Documents/VPA')
setwd('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA') #dir
dir.create('Processed')
ds=read.csv('Info/GSM-Meta_combined_annotated_withGPLs_adjusted_arrays_only.csv')
ds2=ds[grep('Affy',ds$Description),]

gses=unique(ds2$GSE)

list.files('rawData',recursive = T,full.names = T)->files
files=files[grep(paste(gses,collapse = '|'),files)]
files=files[grep('[.]tar',files)]

# there is a difference in the length of 'files' and 'gses; vectors,
# possibly some GSEs are not available, need to check which
fls=do.call(rbind,strsplit(files,'[/]'))[,2]
setdiff(gses,fls) # looks like GSE666 does not have any raw data files

fails=lapply(files,function(x){
  alldr=strsplit(x,'[/]')[[1]]
  exdr=paste(alldr[1:2],collapse='/')
  gse=alldr[2]
  untar(x,exdir = exdr)
  setwd(exdr)
  fls=list.files('.',pattern = 'CEL|cel')
  
  gpls=unique(ds[ds$GSE==gse,]$Description)
  
  fl=lapply(gpls,function(y){
    sub_gsms=ds[ds$Description==y & ds$GSE==gse,]$GSM
    fls2=fls[grep(paste(sub_gsms,collapse='|'),fls)]
    celfiles=try(read.celfiles(fls2))
    if(class(celfiles)=='try-error'){
      fail=c(gse,y,'FAILED')
    }else{
      celfilesNorm=rma(celfiles,normalize=F)
      gplID=strsplit(y,' ')[[1]][1]
      write.exprs(celfilesNorm,paste0('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/Processed/',gse,'_',gplID,'.csv'))
      #write.exprs(celfilesNorm,paste0('~/Documents/VPA/Processed/',gse,'_',gplID,'.csv'))
      fail=c(gse,y,'OK')
    }
    fail
    })
  #setwd('~/Documents/VPA')
  setwd('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA') #dir
  do.call(rbind,fl)
})

fails=do.call(rbind,fails)  
write.csv(fails,'Processed/QC.csv',row.names = F)


# reading GSE5258 that failed using pd.hg.u133a annotations

#meta=read.csv('~/Documents/VPA/Info/GSM-Meta_combined_annotated_withGPLs_adjusted_arrays_only.csv')
meta=read.csv('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/Info/GSM-Meta_combined_annotated_withGPLs_adjusted_arrays_only.csv')
setwd('/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/rawData/GSE5258')
gpl3921=meta[meta$GSE=='GSE5258' & meta$GPL=="'GPL3921'",]$GSM
gpl3921=paste(gpl3921,collapse = '|')
print(gpl3921)

fl=list.files('.',pattern = 'CEL')

celfiles=read.celfiles(fl[grep(gpl3921,fl)],pkgname = 'pd.hg.u133a')

celfilesNorm=rma(celfiles,normalize=F)
#write.exprs(celfilesNorm,'~/Documents/VPA/Processed/GSE5258_[HT_HG-U133A].csv')
write.exprs(celfilesNorm,'/home/jaic1/Documents/GitRepo_Practiceals/IRP/VPA/Processed/GSE5258_[HT_HG-U133A].csv')