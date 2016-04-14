library(AnalyzeFMRI)
library(ggplot2)
library(spatstat)

folder <-  "~/Documents/Onderzoek/Studie_meta_analysis/EffectSizeResults/"

setwd(folder)

results_file <-paste(folder,"peagGTn20.txt",sep="")
res10_file <- paste(folder,"resultsn10.txt",sep="")

results <- read.table(results_file,header=TRUE)
res10 <- read.table(res10_file,header=TRUE)
res10$Cluster <- factor(res10$Cluster)

ref <- f.read.nifti.volume("ClustIt1.nii")[,,,1]
reftab <- melt(ref)
names(reftab) <- c("x","y","z","value")
reftab$bin <- ifelse(reftab$value<4,reftab$value,0)
reftab$bin <- ifelse(reftab$bin==0,0,1)

for(z in min(reftab$z):max(reftab$z)){cat(z,":",sum(reftab$bin[reftab$z==z]==1),"\n")}

MNI <- f.read.nifti.volume("MNI152_T1_2mm_brain.nii")[,,,1]
MNItab <- melt(MNI)
names(MNItab) <- c("x","y","z","value")


pdf("peaksfig.pdf",height=14,width=11)
for(zsel in c(27,40,50)){
MNIsel = MNItab[MNItab$z==zsel,]
subset <- res10[res10$zcord==zsel,]
subcl <- reftab[reftab$z==zsel,]
subcl <- subcl[subcl$bin==1,]

p <- ggplot() + geom_tile(data=MNIsel,aes(x=x,y=y,fill=value)) + scale_fill_gradient(low="grey70",high="grey90")
p <- p + geom_tile(data=subcl,aes(x=x,y=y))
p <- p + geom_point(data=subset,aes(x=xcord,y=ycord,colour=Peakt))+ scale_colour_gradient(low = "red",high="yellow",limits=c(5,20)) + geom_jitter()
print(p)
}
dev.off()