library(oro.nifti)

row.match<-function (x, table, nomatch = NA)
{
    if (class(table) == "matrix")
        table <- as.data.frame(table)
    if (is.null(dim(x)))
        x <- as.data.frame(matrix(x, nrow = 1))
    cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
    ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
    match(cx, ct, nomatch = nomatch)
}

#####

args = commandArgs(trailingOnly=TRUE)

HPCfolder <- "/user/data/gent/gvo000/gvo00022/00_data/05_motor/"
scratchdir <- "/user/scratch/gent/vsc418/vsc41855/Tempresultsn10/"
resultsdir<-"/user/scratch/gent/vsc418/vsc41855/resultsn10/"

#####

copefile <- paste(HPCfolder,"MOTOR_13_cope.nii.gz",sep="")
cope1 <- readNIfTI(copefile)[,,,1:181]

######

setwd(scratchdir)

nfull<-101
nclust<-4
nsample<-10
nsim<-3
nstudy<-4
ngroundtruth<-c(1:181)

it<-as.numeric(args[1])
set.seed(it)

permn1<-sample(ngroundtruth,length(ngroundtruth),replace=FALSE)
gt<-permn1[1:nfull]
permn<-permn1[(nfull+1):(nfull+nstudy*nsample)]

study1<-permn[1:nsample]
study2<-permn[(nsample+1):(2*nsample)]
study3<-permn[(2*nsample+1):(3*nsample)]
study4<-permn[(3*nsample+1):(4*nsample)]
studies<-cbind(study1,study2,study3,study4)

cope<-cope1[,,,gt]

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
pmap <- 1-pt(tmap,nfull-1)

pname <- paste(scratchdir,"pstat1",sep="")
writeNIfTI(pmap,filename=pname)

zname <- paste(scratchdir,"zmap",sep="")
zmap<-qnorm(pt(tmap,nfull-1))
writeNIfTI(zmap,zname)

rescope<- array(NA,dim=c(91,109,91,nfull))
avcope<-apply(cope,c(1,2,3),mean)
for(i in 1:nfull)
{rescope[,,,i]<-cope[,,,i]-avcope}

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
imap2<-length(clusterresults$index)-imap+1
SaveClustGT<-paste(resultsdir,"ClustIt",it,sep="")
writeNIFTI(imap2,SaveClustGT)

xco<-read.table("peakfile.txt",skip=1)$V3
yco<-read.table("peakfile.txt",skip=1)$V4
zco<-read.table("peakfile.txt",skip=1)$V5

peakresults<-read.table("peakfile.txt",skip=1)
peakresults2<-cbind(rep(it,dim(peakresults)[1]),peakresults,tmap[cbind(xco+1,yco+1,zco+1)])
write.table(peakresults2,paste(resultsdir,"peakGT.txt",sep=""),append=T,col.names=FALSE,row.names=FALSE)

for(st in 1:nstudy)
{
print(st)

copes1<-cope1[,,,studies[,st]]

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

clusterindexStudy<-read.table(clusterfile1,skip=1)


indicesStudy<-paste(scratchdir,"indexfiles1.nii.gz",sep="")
imapStudy<-readNIfTI(indicesStudy)[,,]


for(cl in 1:nclust)
{
  maskc1<-imap==clusterresults$index[cl]
compc1<-which(peakmasks1==1&maskc1==1,arr.ind=T)-1
cord<-data.frame(cbind(xcord,ycord,zcord))
a1<-row.match(data.frame(compc1),cord)
boolean<-ifelse(length(a1)!=0,1,0)
peakmap<-readNIfTI(paste(scratchdir,"outputs",st,"_tstat1.nii.gz",sep=""))[,,]
es<-mean(tmap[maskc1==1]/sqrt(nfull))
wes<-es*(1-3/(4*(nfull-1)-1))
if(boolean==1)
{
clind<-read.table("peakfiles1.txt",skip=1)$V1[a1]
maskSt<-array(0,dim=c(91,109,91))
for(tel in 1:length(unique(clind)))
{maskSt[imapStudy==clind[tel]]<-1}
doorsnede<-which(maskSt==1&maskc1==1,arr.ind=T)
pes<-mean(tmap[doorsnede]/sqrt(nfull))
peakt<-peakmap[compc1+1]
peak<-peakmap[compc1+1]/sqrt(nsample)
wpeak<-peak*(1-3/(4*(nsample-1)-1))
tabres<-cbind(rep(it,length(a1)),rep(st,length(a1)),peakt,peak,wpeak,rep(cl,length(a1)),rep(boolean,length(a1)),rep(es,length(a1)),rep(wes,length(a1)),rep(pes,length(a1)),xcord[a1],ycord[a1],zcord[a1])
}
if(boolean==0)
{
peak<-max(peakmap[maskc1==1])/sqrt(nsample)
compmax<-which(peakmap==max(peakmap[maskc1==1]),arr.ind=T)
pes<-es
wpeak<-peak*(1-3/(4*(nsample-1)-1))
tabres<-cbind(it,st,max(peakmap[maskc1==1]),peak,wpeak,cl,boolean,es,wes,pes,compmax[1],compmax[2],compmax[3])
}
write.table(tabres,paste(resultsdir,"results.txt",sep=""),append=T,col.names=FALSE,row.names=FALSE)
}


}
