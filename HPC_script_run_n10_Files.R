HPCfolder <- "/user/data/gent/gvo000/gvo00022/00_data/05_motor/"
scratchdir <- "/user/scratch/gent/vsc418/vsc41855/Tempresultsn10/"
resultsdir<-"/user/scratch/gent/vsc418/vsc41855/resultsn10/"


setwd(scratchdir)

write(c("Iteration","Study","Peakt","Peak","WPeak","Cluster","Significant","ESC","WESC","PES","xcord","ycord","zcord"),paste(resultsdir,"results.txt",sep=""),ncolumns=12)
write(c("Iteration","Cluster","Peakz","xcord","ycord","zcord","Peakt"),paste(resultsdir,"peakGT.txt",sep=""),ncolumns=6)
