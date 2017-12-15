args <- commandArgs(TRUE)
sDGEpath<-args[1]
mDGEpath<-args[2]
Samplename<-args[3]

# library(scater)

setwd(sDGEpath)
filelist<-list.files("./",pattern=".txt")
formfile<-read.table(filelist[1],header=T,row.names=1)
#init a data frame
DGE<-data.frame(formfile[6])
makeDge<-function(name){
	xxx<-read.table(name,header=T,row.names=1)
	value<-as.data.frame(xxx[6])
	names(value)<-paste(unlist(strsplit(name,"_"))[1:2],collapse="_")
	return(value)
}
DGE<-makeDge(filelist[1])
for(i in 2:length(filelist)){
	new<-makeDge(filelist[i])
	DGE<-data.frame(DGE,new)
}
outputname<-paste(mDGEpath,"/",Samplename,".Rdata",sep="")
save(DGE,file=outputname)