library(oro.nifti)
library(prodlim)

HOMEfolder <- "/Users/bmoerker/Documents/server/fMRI/Stanford/Softwarecode/"
HPCfolder <- "/Users/bmoerker/Desktop/contrast/80_Unrelated_cope/"
scratchdir <- "/Users/bmoerker/Documents/server/fMRI/Stanford/Softwarecode/Scratch/"
resultsdir<-"/Users/bmoerker/Documents/server/fMRI/Stanford/Softwarecode/Results/"

setwd(scratchdir)

write(c("Iteration","Study","Peak","WPeak","Cluster","Significant","ESC","WESC"),paste(resultsdir,"results.txt",sep=""),ncolumns=8)

## Reading in "full data" to create ground truth
## For now, 30 subjects as ground truth - sample from other 50 subjects

#####################################################################
## Identifying clusters using super-ground truth

nindiv<-78

cope <- array(NA,dim=c(91,109,91,nindiv))
for (s in 1:nindiv){
copefile <- paste(HPCfolder,"subject_",s,"_contrast_1",sep="")
cope1 <- readNIfTI(copefile)[,,]
cope[,,,s] <- cope1
}

mask_subject <- ifelse(cope==0,0,1)
mask_av <- apply(mask_subject,c(1,2,3),mean)
mask <- ifelse(mask_av==1,1,0)
maskname <- paste(scratchdir,"mask",sep="")
writeNIfTI(mask,maskname)

copename <- paste(scratchdir,"cope",sep="")
writeNIfTI(cope,filename=copename)

fsl_model_cmd = paste("randomise -i ",copename," -o output -1 --glm_output",sep="")
system(fsl_model_cmd)

tname <- paste(scratchdir,"output_tstat1.nii.gz",sep="")
tmap <- readNIfTI(tname)[,,]
pmap <- 1-pt(tmap,nindiv-1)

pname <- paste(scratchdir,"pstat1",sep="")
writeNIfTI(pmap,filename=pname)

zname <- paste(scratchdir,"zmap",sep="")
zmap<-qnorm(pt(tmap,nindiv-1))
writeNIfTI(zmap,zname)

rescope<- array(NA,dim=c(91,109,91,nindiv))
avcope<-apply(cope,c(1,2,3),mean)
for(i in 1:nindiv)
{rescope[,,,i]<-cope[,,,i]-avcope}
  resname <- paste(scratchdir,"rescope",sep="")
  writeNIfTI(rescope,filename=resname)

df<-nindiv-1
smoothfile <- paste(scratchdir,"smoothness.txt",sep="")
system(paste("smoothest -d ",df," -m ",maskname," -r ",resname, ">",smoothfile,sep=""))
smooth <- read.table(smoothfile)
vol <- smooth$V2[2]
dlh <- smooth$V2[1]
res <- smooth$V2[3]

RESELcount<-vol/res
FWE_cmd<-paste("ptoz 0.05 -g ",RESELcount)
FWEthresh<-as.numeric(system(FWE_cmd,intern=TRUE))

sign <- ifelse((zmap>FWEthresh),1,0)

clusterfile <- paste(scratchdir,"clusterfile.txt",sep="")
cl_cmd <- paste("cluster --in=",zname," --thresh=",FWEthresh," -p ",1," -d ",dlh," -o cluster --volume=",vol," --olmax=peakfile.txt --oindex=indexfile > ",clusterfile,sep="")
system(cl_cmd)

clusterresults <- read.table(clusterfile,skip=1)
names(clusterresults) <- c("index","voxels","p","logp","max","x","y","z","cogx","cogy","cogz")

indices<-paste(scratchdir,"indexfile.nii.gz",sep="")

imapSGT<-readNIfTI(indices)[,,]
nclustSGT<-length(clusterresults$index)
IndexSGT0<-clusterresults$index
IndexSGT<-rev(clusterresults$index)


#################




nfull<-30
nclust<-3
nsample<-12
nsim<-3
nstudy<-4
ngroundtruth<-c(1:78)

for(it in 1:nsim)
{
print(it)

set.seed(it)
permn1<-sample(ngroundtruth,length(ngroundtruth),replace=FALSE)
gt<-permn1[1:nfull]
permn<-permn1[(nfull+1):(nfull+nstudy*nsample)]

study1<-permn[1:nsample]
study2<-permn[(nsample+1):(2*nsample)]
study3<-permn[(2*nsample+1):(3*nsample)]
study4<-permn[(3*nsample+1):(4*nsample)]
studies<-cbind(study1,study2,study3,study4)


cope <- array(NA,dim=c(91,109,91,nfull))
for (s in 1:nfull){
copefile <- paste(HPCfolder,"subject_",gt[s],"_contrast_1",sep="")
cope1 <- readNIfTI(copefile)[,,]
cope[,,,s] <- cope1
}

## Create mask before analyses
## ? use also this mask for future analyses? - now rerun each time


mask_subject <- ifelse(cope==0,0,1)
mask_av <- apply(mask_subject,c(1,2,3),mean)
mask <- ifelse(mask_av==1,1,0)
maskname <- paste(scratchdir,"mask",sep="")
writeNIfTI(mask,maskname)

copename <- paste(scratchdir,"cope",sep="")
writeNIfTI(cope,filename=copename)

## Ground truth: FWE analysis with control at level 0.05

fsl_model_cmd = paste("randomise -i ",copename," -o output -1 --glm_output",sep="")
system(fsl_model_cmd)

tname <- paste(scratchdir,"output_tstat1.nii.gz",sep="")
tmap <- readNIfTI(tname)[,,]
pmap <- 1-pt(tmap,nfull-1)

pname <- paste(scratchdir,"pstat1",sep="")
writeNIfTI(pmap,filename=pname)

zname <- paste(scratchdir,"zmap",sep="")
zmap<-qnorm(pt(tmap,nfull-1))
writeNIfTI(zmap,zname)

# residuals for smoothing

rescope<- array(NA,dim=c(91,109,91,nfull))
avcope<-apply(cope,c(1,2,3),mean)
for(i in 1:nfull)
{
  rescope[,,,i]<-cope[,,,i]-avcope}
  resname <- paste(scratchdir,"rescope",sep="")
  writeNIfTI(rescope,filename=resname)

df<-nfull-1
smoothfile <- paste(scratchdir,"smoothness.txt",sep="")
system(paste("smoothest -d ",df," -m ",maskname," -r ",resname, ">",smoothfile,sep=""))
smooth <- read.table(smoothfile)
vol <- smooth$V2[2]
dlh <- smooth$V2[1]
res <- smooth$V2[3]


RESELcount<-vol/res
FWE_cmd<-paste("ptoz 0.05 -g ",RESELcount)

FWEthresh<-as.numeric(system(FWE_cmd,intern=TRUE))

# how many voxels are active?
sign <- ifelse((zmap>FWEthresh),1,0)
sum(sign==1)

# How many clusters? -> use cluster command

clusterfile <- paste(scratchdir,"clusterfile.txt",sep="")
cl_cmd <- paste("cluster --in=",zname," --thresh=",FWEthresh," -p ",1," -d ",dlh," -o cluster --volume=",vol," --olmax=peakfile.txt --oindex=indexfile > ",clusterfile,sep="")
system(cl_cmd)

clusterresults <- read.table(clusterfile,skip=1)
names(clusterresults) <- c("index","voxels","p","logp","max","x","y","z","cogx","cogy","cogz")

indices<-paste(scratchdir,"indexfile.nii.gz",sep="")
imap<-readNIfTI(indices)[,,]


clustindicator<-vector("numeric",3)
for(ell in 1:3)
{a2<-vector("numeric",length(nclustSGT))
clust<-clusterresults$index[ell]
clust1<-imap==clust
for(kell in 1:nclustSGT)
{clust2<-imapSGT==IndexSGT0[kell]
 a2[kell]<-sum(clust1&clust2)
  }
a1<-which(a2==max(a2))[1]
clustindicator[ell]<-IndexSGT[a1]
}


# Retain three largest clusters

# clust1<-clusterresults$index[1]
# clust2<-clusterresults$index[2]
# clust3<-clusterresults$index[3]
#
# maskc1<-imap==clust1
# maskc2<-imap==clust2
# maskc3<-imap==clust3
#
# csize1<-sum(maskc1)
# csize2<-sum(maskc2)
# csize3<-sum(maskc3)
#
# # effect sizes: using t or z? with z, infinity values arise
#
# es1<-mean(tmap[maskc1==1]/sqrt(nfull))
# es2<-mean(tmap[maskc2==1]/sqrt(nfull))
# es3<-mean(tmap[maskc3==1]/sqrt(nfull))
#
# # weighted effect sizes, cfr. Cohen's d
#
# J1<-1-3/(4*(nfull-1)-1)
#
# wes1<-es1*J1
# wes2<-es2*J1
# wes3<-es3*J1


#############



for(st in 1:nstudy)
{
print(st)

copes1 <- array(NA,dim=c(91,109,91,nsample))
for (s in 1:nsample){
copefile1 <- paste(HPCfolder,"subject_",studies[s,st],"_contrast_1",sep="")
cope1 <- readNIfTI(copefile1)[,,]
copes1[,,,s] <- cope1
}

copename1 <- paste(scratchdir,"copes1",sep="")
writeNIfTI(copes1,filename=copename1)

rescope1<- array(NA,dim=c(91,109,91,nsample))
avcope1<-apply(copes1,c(1,2,3),mean)
for(i in 1:nsample)
{
rescope1[,,,i]<-copes1[,,,i]-avcope1
}
resname1 <- paste(scratchdir,"rescope1",sep="")
writeNIfTI(rescope1,filename=resname1)

mask_subject <- ifelse(copes1==0,0,1)
mask_av <- apply(mask_subject,c(1,2,3),mean)
mask <- ifelse(mask_av==1,1,0)
maskname1 <- paste(scratchdir,"mask1",sep="")
writeNIfTI(mask,maskname1)

fsl_model_cmd = paste("randomise -i ",copename1," -o outputs",st," -1 --glm_output",sep="")
system(fsl_model_cmd)
tname1 <- paste(scratchdir,"outputs",st,"_tstat1.nii.gz",sep="")
tmap1 <- readNIfTI(tname1)[,,]
zname1 <- paste(scratchdir,"zmap1",sep="")
zmap1<-qnorm(pt(tmap1,nsample-1))
writeNIfTI(zmap1,zname1)

df<-nsample-1
smoothfile1 <- paste(scratchdir,"smoothness1.txt",sep="")
system(paste("smoothest -d ",df," -m ",maskname1," -r ",resname1, ">",smoothfile1,sep=""))
smooth <- read.table(smoothfile1)
vol <- smooth$V2[2]
dlh <- smooth$V2[1]
clusterfile1 <- paste(scratchdir,"clusterfile1.txt",sep="")
cl_cmd <- paste("cluster --in=",zname1," --thresh=",qnorm(1-0.001)," -p ",0.01," -d ",dlh," -o cluster --volume=",vol," --olmax=peakfiles1.txt --oindex=indexfiles1 > ",clusterfile1,sep="")
system(cl_cmd)


xcord<-read.table("peakfiles1.txt",skip=1)$V3
ycord<-read.table("peakfiles1.txt",skip=1)$V4
zcord<-read.table("peakfiles1.txt",skip=1)$V5

peakmasks1 <- array(0,dim=c(91,109,91))
peakmasks1[cbind(xcord+1,ycord+1,zcord+1)]<-1


for(cl in 1:nclust)
{maskc1<-imap==clusterresults$index[cl]
compc1<-which(peakmasks1==1&maskc1==1,arr.ind=T)-1
cord<-data.frame(cbind(xcord,ycord,zcord))
a1<-row.match(data.frame(compc1),cord)
boolean<-ifelse(length(a1)!=0,1,0)
peakmap<-readNIfTI(paste(scratchdir,"outputs",st,"_tstat1.nii.gz",sep=""))[,,]
es<-mean(tmap[maskc1==1]/sqrt(nfull))
wes<-es*(1-3/(4*(nfull-1)-1))
if(boolean==1)
{
peak<-peakmap[compc1+1]/sqrt(nsample)
wpeak<-peak*(1-3/(4*(nsample-1)-1))
tabres<-cbind(rep(it,length(a1)),rep(st,length(a1)),peak,wpeak,rep(cl,length(a1)),rep(boolean,length(a1)),rep(es,length(a1)),rep(wes,length(a1)))
}
if(boolean==0)
{
peak<-max(peakmap[maskc1==1])/sqrt(nsample)
wpeak<-peak*(1-3/(4*(nsample-1)-1))
tabres<-cbind(it,st,peak,wpeak,cl,boolean,es,wes)
}
write.table(tabres,paste(resultsdir,"results.txt",sep=""),append=T,col.names=FALSE,row.names=FALSE)
}


}


}
