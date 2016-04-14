library(AnalyzeFMRI)
library(ggplot2)

folder <-  "~/Documents/Onderzoek/Studie_meta_analysis/EffectSizeResults/"

setwd(folder)

results_file <-paste(folder,"peagGTn20.txt",sep="")

results <- read.table(results_file,header=TRUE)

#MNI <- f.read.nifti.volume("MNI152_T1_2mm_brain.nii")[,,,1]

pl <- ggplot(results,aes(x=xcord,y=ycord,colour=Peakt,colour=Peakt))

pl + geom_point() + geom_jitter()


