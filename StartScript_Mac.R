library(oro.nifti)

HOMEfolder <- "/Users/bmoerker/Documents/server/fMRI/Stanford/Softwarecode/"
HPCfolder <- "/Users/bmoerker/Desktop/contrast/80_Unrelated_cope/"
scratchdir <- "/Users/bmoerker/Documents/server/fMRI/Stanford/Softwarecode/Scratch/"

setwd(scratchdir)

## Reading in "full data" to create ground truth
## For now, 30 subjects as ground truth - sample from other 50 subjects

nfull <- 30
cope <- array(NA,dim=c(91,109,91,30))
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

# smoothness op z-statistieken? -> geeft extreme estimates

smoothfile <- paste(scratchdir,"smoothness.txt",sep="")
system(paste("smoothest -z ",tname," -m ",maskname,">",smoothfile,sep=""))
smooth <- read.table(smoothfile)
vol <- smooth$V2[2]
dlh <- smooth$V2[1]
res <- smooth$V2[3]
clusterfile <- paste(scratchdir,"clusterfile.txt",sep="")

RESELcount<-vol/res
FWE_cmd<-paste("ptoz 0.05 -g ",RESELcount)

FWEthresh<-as.numeric(system(FWE_cmd,intern=TRUE))

# Code Joke -> komt niet helemaal overeen
#Res<-RESELcount
#Res<-sum(mask)/(5/2)**3
#zs<-seq(1,15,0.0001)
#pN_RFT<-Res*exp(-zs**2/2)*zs**2
#cutoff_RFT<-min(zs[pN_RFT<0.05])


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


#############
