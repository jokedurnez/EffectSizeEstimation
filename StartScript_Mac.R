library(oro.nifti)
library(prodlim)

HOMEfolder <- "/Users/bmoerker/Documents/server/fMRI/Stanford/Softwarecode/"
HPCfolder <- "/Users/bmoerker/Desktop/contrast/80_Unrelated_cope/"
scratchdir <- "/Users/bmoerker/Documents/server/fMRI/Stanford/Softwarecode/Scratch/"

setwd(scratchdir)

write(c("Study","Peak","Cluster","Significant","ESC"),paste(scratchdir,"results.txt",sep=""),ncolumns=5)

## Reading in "full data" to create ground truth
## For now, 30 subjects as ground truth - sample from other 50 subjects

nfull <- 30
cope <- array(NA,dim=c(91,109,91,nfull))
for (s in 1:nfull){
copefile <- paste(HPCfolder,"subject_",s,"_contrast_1",sep="")
cope1 <- readNIfTI(copefile)[,,]
cope[,,,s] <- cope1
}

## Create mask before analyses
## ? use also this mask for future analyses?
## Is only used for command cluster, not for ground truth?

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
{ print(i)
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

# Retain three largest clusters

clust1<-clusterresults$index[1]
clust2<-clusterresults$index[2]
clust3<-clusterresults$index[3]

maskc1<-imap==clust1
maskc2<-imap==clust2
maskc3<-imap==clust3

csize1<-sum(maskc1)
csize2<-sum(maskc2)
csize3<-sum(maskc3)

# effect sizes: using t or z? with z, infinity values arise

es1<-mean(tmap[maskc1==1]/sqrt(nfull))
es2<-mean(tmap[maskc2==1]/sqrt(nfull))
es3<-mean(tmap[maskc3==1]/sqrt(nfull))

# weighted effect sizes, cfr. Cohen's d

J1<-1-3/(4*(nfull-1)-1)

wes1<-es1*J1
wes2<-es2*J1
wes3<-es3*J1


#############

ngroundtruth<-c(31:78)

permn<-sample(ngroundtruth,length(ngroundtruth),replace=FALSE)
nsample<-12
nstudy<-4
study1<-permn[1:nsample]
study2<-permn[(nsample+1):(2*nsample)]
study3<-permn[(2*nsample+1):(3*nsample)]
study4<-permn[(3*nsample+1):(4*nsample)]


copes1 <- array(NA,dim=c(91,109,91,nsample))
copes2 <- array(NA,dim=c(91,109,91,nsample))
copes3 <- array(NA,dim=c(91,109,91,nsample))
copes4 <- array(NA,dim=c(91,109,91,nsample))
for (s in 1:nsample){
copefile1 <- paste(HPCfolder,"subject_",study1[s],"_contrast_1",sep="")
copefile2 <- paste(HPCfolder,"subject_",study2[s],"_contrast_1",sep="")
copefile3 <- paste(HPCfolder,"subject_",study3[s],"_contrast_1",sep="")
copefile4 <- paste(HPCfolder,"subject_",study4[s],"_contrast_1",sep="")
cope1 <- readNIfTI(copefile1)[,,]
cope2 <- readNIfTI(copefile2)[,,]
cope3 <- readNIfTI(copefile3)[,,]
cope4 <- readNIfTI(copefile4)[,,]
copes1[,,,s] <- cope1
copes2[,,,s] <- cope2
copes3[,,,s] <- cope3
copes4[,,,s] <- cope4
}

copename1 <- paste(scratchdir,"copes1",sep="")
writeNIfTI(copes1,filename=copename)
copename2 <- paste(scratchdir,"copes2",sep="")
writeNIfTI(copes2,filename=copename)
copename3 <- paste(scratchdir,"copes3",sep="")
writeNIfTI(copes3,filename=copename)
copename4 <- paste(scratchdir,"copes4",sep="")
writeNIfTI(copes4,filename=copename)

rescope1<- array(NA,dim=c(91,109,91,nsample))
rescope2<- array(NA,dim=c(91,109,91,nsample))
rescope3<- array(NA,dim=c(91,109,91,nsample))
rescope4<- array(NA,dim=c(91,109,91,nsample))
avcope1<-apply(copes1,c(1,2,3),mean)
avcope2<-apply(copes2,c(1,2,3),mean)
avcope3<-apply(copes3,c(1,2,3),mean)
avcope4<-apply(copes4,c(1,2,3),mean)
for(i in 1:nsample)
{ print(i)
  rescope1[,,,i]<-copes1[,,,i]-avcope1
rescope2[,,,i]<-copes2[,,,i]-avcope2
rescope3[,,,i]<-copes3[,,,i]-avcope3
rescope4[,,,i]<-copes4[,,,i]-avcope4}
resname1 <- paste(scratchdir,"rescope1",sep="")
writeNIfTI(rescope1,filename=resname1)
resname2 <- paste(scratchdir,"rescope2",sep="")
writeNIfTI(rescope2,filename=resname2)
resname3 <- paste(scratchdir,"rescope3",sep="")
writeNIfTI(rescope3,filename=resname3)
resname4 <- paste(scratchdir,"rescope4",sep="")
writeNIfTI(rescope4,filename=resname4)

mask_subject <- ifelse(copes1==0,0,1)
mask_av <- apply(mask_subject,c(1,2,3),mean)
mask <- ifelse(mask_av==1,1,0)
maskname1 <- paste(scratchdir,"mask1",sep="")
writeNIfTI(mask,maskname1)

mask_subject <- ifelse(copes2==0,0,1)
mask_av <- apply(mask_subject,c(1,2,3),mean)
mask <- ifelse(mask_av==1,1,0)
maskname2 <- paste(scratchdir,"mask2",sep="")
writeNIfTI(mask,maskname2)

mask_subject <- ifelse(copes3==0,0,1)
mask_av <- apply(mask_subject,c(1,2,3),mean)
mask <- ifelse(mask_av==1,1,0)
maskname3 <- paste(scratchdir,"mask3",sep="")
writeNIfTI(mask,maskname3)

mask_subject <- ifelse(copes4==0,0,1)
mask_av <- apply(mask_subject,c(1,2,3),mean)
mask <- ifelse(mask_av==1,1,0)
maskname4 <- paste(scratchdir,"mask4",sep="")
writeNIfTI(mask,maskname4)

fsl_model_cmd = paste("randomise -i ",copename1," -o outputs1 -1 --glm_output",sep="")
system(fsl_model_cmd)
tname1 <- paste(scratchdir,"outputs1_tstat1.nii.gz",sep="")
tmap1 <- readNIfTI(tname1)[,,]
zname1 <- paste(scratchdir,"zmap1",sep="")
zmap1<-qnorm(pt(tmap1,nsample-1))
writeNIfTI(zmap1,zname1)

fsl_model_cmd = paste("randomise -i ",copename2," -o outputs2 -1 --glm_output",sep="")
system(fsl_model_cmd)
tname2 <- paste(scratchdir,"outputs2_tstat1.nii.gz",sep="")
tmap2 <- readNIfTI(tname2)[,,]
zname2 <- paste(scratchdir,"zmap2",sep="")
zmap2<-qnorm(pt(tmap2,nsample-1))
writeNIfTI(zmap2,zname2)

fsl_model_cmd = paste("randomise -i ",copename3," -o outputs3 -1 --glm_output",sep="")
system(fsl_model_cmd)
tname3 <- paste(scratchdir,"outputs3_tstat1.nii.gz",sep="")
tmap3 <- readNIfTI(tname3)[,,]
zname3 <- paste(scratchdir,"zmap3",sep="")
zmap3<-qnorm(pt(tmap3,nsample-1))
writeNIfTI(zmap3,zname3)


fsl_model_cmd = paste("randomise -i ",copename4," -o outputs4 -1 --glm_output",sep="")
system(fsl_model_cmd)
tname4 <- paste(scratchdir,"outputs4_tstat1.nii.gz",sep="")
tmap4 <- readNIfTI(tname4)[,,]
zname4 <- paste(scratchdir,"zmap4",sep="")
zmap4<-qnorm(pt(tmap4,nsample-1))
writeNIfTI(zmap4,zname4)

df<-nsample-1
smoothfile1 <- paste(scratchdir,"smoothness1.txt",sep="")
smoothfile2 <- paste(scratchdir,"smoothness2.txt",sep="")
smoothfile3 <- paste(scratchdir,"smoothness3.txt",sep="")
smoothfile4 <- paste(scratchdir,"smoothness4.txt",sep="")

nstudies<-4
nclust<-3



system(paste("smoothest -d ",df," -m ",maskname1," -r ",resname1, ">",smoothfile1,sep=""))
smooth <- read.table(smoothfile1)
vol <- smooth$V2[2]
dlh <- smooth$V2[1]
clusterfile1 <- paste(scratchdir,"clusterfile1.txt",sep="")
cl_cmd <- paste("cluster --in=",zname1," --thresh=",qnorm(1-0.001)," -p ",0.05," -d ",dlh," -o cluster --volume=",vol," --olmax=peakfiles1.txt --oindex=indexfiles1 > ",clusterfile1,sep="")
system(cl_cmd)

xcord<-read.table("peakfiles1.txt",skip=1)$V3
ycord<-read.table("peakfiles1.txt",skip=1)$V4
zcord<-read.table("peakfiles1.txt",skip=1)$V5


studynr<-1



clusternr<-1

peakmasks1 <- array(0,dim=c(91,109,91))
peakmasks1[cbind(xcord+1,ycord+1,zcord+1)]<-1
st<-1

for(cl in 1:nclust)
{print(cl)
  maskc1<-imap==clusterresults$index[cl]
compc1<-which(peakmasks1==1&maskc1==1,arr.ind=T)-1
cord<-data.frame(cbind(xcord,ycord,zcord))
a1<-row.match(data.frame(compc1),data.frame(cbind(xcord,ycord,zcord)))
boolean<-ifelse(length(a1)!=0,1,0)
peakmap<-readNIfTI(paste(scratchdir,"output_tstat",st,".nii.gz",sep=""))[,,]
es<-mean(tmap[maskc1==1]/sqrt(nfull))
if(boolean==1)
{

peak<-peakmap[compc1+1]/sqrt(nsample)
tabres<-cbind(rep(studynr,length(a1)),peak,rep(cl,length(a1)),rep(boolean,length(a1)),rep(es,length(a1)))
}
if(boolean==0)
{
peak<-max(peakmap[maskc1==1])/sqrt(nsample)
tabres<-cbind(studynr,peak,cl,boolean,es)
}
write(t(tabres),paste(scratchdir,"results.txt",sep=""),append=T)
}
