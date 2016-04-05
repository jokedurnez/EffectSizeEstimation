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

## Ground truth: FDR analysis with q=0.05

fsl_model_cmd = paste("randomise -i ",copename," -o output -1 --glm_output",sep="")
system(fsl_model_cmd)

tname <- paste(scratchdir,"output_tstat1.nii.gz",sep="")
tmap <- readNIfTI(tname)[,,]
pmap <- 1-pt(tmap,nfull-1)

pname <- paste(scratchdir,"pstat1",sep="")
writeNIfTI(pmap,filename=pname)

thresfile <- paste(scratchdir,"threshold.txt",sep="")
fdr_cmd <- paste("fdr -i",pname,"-q 0.05 >",thresfile)
system(fdr_cmd)

fdr_threshold <- read.table(thresfile,skip=1)$V1

# how many voxels are active?
sign <- ifelse(pmap<fdr_threshold,1,0)
sum(sign==1)

# How many clusters? -> use cluster command

# First: transform from t to z

zmap<-qnorm(pt(tmap,nfull-1))

fdr_threshold_t<-qt(1-fdr_threshold,nfull-1)
fdr_threshold_z<-qnorm(pt(fdr_threshold_t,nfull-1))

zname <- paste(scratchdir,"zmap",sep="")
writeNIfTI(zmap,zname)

smoothfile <- paste(scratchdir,"smoothness.txt",sep="")
system(paste("smoothest -z ",tname," -m ",maskname,">",smoothfile,sep=""))
smooth <- read.table(smoothfile)
vol <- smooth$V2[2]
dlh <- smooth$V2[1]
clusterfile <- paste(scratchdir,"clusterfile.txt",sep="")
cl_cmd <- paste("cluster --in=",zname," --thresh=",fdr_threshold_z," -p ",1," -d ",dlh," -o cluster --volume=",vol," --olmax=peakfile.txt > ",clusterfile,sep="")
system(cl_cmd)

clusterresults <- read.table(clusterfile,skip=1)
names(clusterresults) <- c("index","voxels","p","logp","max","x","y","z","cogx","cogy","cogz")
