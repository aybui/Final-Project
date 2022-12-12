args <- commandArgs(trailingOnly = TRUE)
project <- args[1]
nTaxa <- as.numeric(args[2])
NLOCI <- as.numeric(args[3])
MissingAllowed <- as.numeric(args[4])
folderTail <- args[5]

print("NOTE: PDFs GENERATED HERE ARE BEST VIEWED IN SAFARI OR SIMILAR (NOT IN PREVIEW)")
print("STEPS:")
print(paste(c("(1) Observe ../Results/",project,"_AlignmentSummary_MaxTally",folderTail,".pdf"),sep="",collapse=""))
print("(2) Choose a value for propToBeGood")
print("(3) Rerun  TrimAndMaskRawAlignments")
print("(4) Run PlotAlignmentSummary.r")
print(paste(c("(5) Observe ../Results/",project,"_AlignmentSummary_WindowTally",folderTail,".pdf"),sep="",collapse=""))
print("(6) Choose a value for minGoodSites") 
print("(7) Run TrimAndMaskRawAlignments")
print("(8) Run PlotAlignmentSummary.r")
print("(9) Check trimmed alignments and repeat 1-9 if necessary.")
print("(10) Finally rerun TrimAndMaskRawAlignments and PlotAlignmentSummary.r one final time.")

blah<<-as.matrix(read.table(paste(c("../Results/",project,"_maxTally",folderTail,".txt"),sep="",collapse="")))
pdf(paste(c("../Results/",project,"_AlignmentSummary_MaxTally",folderTail,".pdf"),sep="",collapse=""))
plot((1:ncol(blah))/ncol(blah),log10(apply(blah,2,sum)),type="l",xlab="Occurrence of Most Common Base Within a Site", ylab="Number of Sites",main=project)
dev.off()

blah<<-as.matrix(read.table(paste(c("../Results/",project,"_windowTally",folderTail,".txt"),sep="",collapse="")))
pdf(paste(c("../Results/",project,"_AlignmentSummary_WindowTally",folderTail,".pdf"),sep="",collapse=""))
plot(log10(apply(blah,1,sum)),type="l",xlab="Number of Good Sites in Window.", ylab="Frequency",main=project)
dev.off()

print("NOTE: IF YOU ENCOUNTER A BUG ON THE NEXT STEP, THEN YOU ARE PROBABLY TRYING TO PERFORM TRIMMING SUMMARY ON A NON-ORTHOLOGY ALIGNMENT,WHICH WON'T WORK BECAUSE THERE ARE DIFFERNT NUMBERS OF TAXA IN THE ALIGNMENTS")

blah<<-as.matrix(read.table(paste(c("../Results/",project,"_trimmedLens",folderTail,".txt"),sep="",collapse="")))
pdf(paste(c("../Results/",project,"_AlignmentSummary_Lengths",folderTail,".pdf"),sep="",collapse=""))
plot(-100,-100,xlim=c(0,(nTaxa-1)),ylim=c(0,max(blah)),xlab="Number of Missing Characters Allowed Per Site", ylab="Alignment Length",main=project)
for(i in 1:nTaxa){points(rep(i-1,NLOCI)+runif(NLOCI,-0.25,0.25),blah[,i],cex=0.1)}
points(rep(MissingAllowed-1,NLOCI)+runif(NLOCI,-0.25,0.25),blah[,MissingAllowed],cex=0.1,col="blue")
#arrows(MissingAllowed-1,max(blah),MissingAllowed-1,max(blah-100),angle=30,length=0.05)
dev.off()

blahT<<-as.matrix(read.table(paste(c("../Results/",project,"_charsByTaxon",folderTail,".txt"),sep="",collapse="")))
pdf(paste(c("../Results/",project,"_AlignmentSummary_CharsByTaxon",folderTail,".pdf"),sep="",collapse=""))
plot(-100,-100,xlim=c(0,(nTaxa-1)),ylim=c(0,max(blahT)),xlab="Number of Missing Characters Allowed Per Site", ylab="Total Number of Non-Missing Characters",main=project)
for(i in 1:nTaxa){points(rep(i-1,nTaxa)+runif(nTaxa,-0.25,0.25),blahT[,i],cex=0.1)}
points(rep(MissingAllowed-1,nTaxa)+runif(nTaxa,-0.25,0.25),blahT[,MissingAllowed],cex=0.1,col="blue")
#arrows(MissingAllowed-1,max(blahT),MissingAllowed-1,max(blahT-10),angle=30,length=0.05)
text(rep(nTaxa+1,nTaxa),blahT[,nTaxa],1:nTaxa,cex=0.5)
dev.off()


pdf(paste(c("../Results/",project,"_AlignmentSummary_nSites",folderTail,".pdf"),sep="",collapse=""))
plot(-100,-100,xlim=c(0,(nTaxa-1)),ylim=c(0,max(apply(blah,2,sum))),xlab="Number of Missing Characters Allowed Per Site", ylab="Total Number of Sites",main=project)
points(0:(nTaxa-1),apply(blah,2,sum),pch=19,cex=0.5)
points(MissingAllowed-1,apply(blah,2,sum)[MissingAllowed],pch=19,cex=1,col="blue")
text(MissingAllowed+5,apply(blah,2,sum)[MissingAllowed],apply(blah,2,sum)[MissingAllowed],cex=0.5,col="blue")
dev.off()

pdf(paste(c("../Results/",project,"_AlignmentSummary_nLoci",folderTail,".pdf"),sep="",collapse=""))
plot(-100,-100,xlim=c(-6,(nTaxa-1)),ylim=c(0,NLOCI),xlab="Number of Missing Characters Allowed Per Site", ylab="Number of Loci",main=project)
lines(c(MissingAllowed-1,MissingAllowed-1),c(0,NLOCI),lwd=5,col="gray")
nLoci=matrix(0,nrow=5,ncol=nTaxa)
for(i in 1:nTaxa){nLoci[1,i]=sum(blah[,i]>1500)}
for(i in 1:nTaxa){nLoci[2,i]=sum(blah[,i]>1000)}
for(i in 1:nTaxa){nLoci[3,i]=sum(blah[,i]>500)}
for(i in 1:nTaxa){nLoci[4,i]=sum(blah[,i]>250)}
for(i in 1:nTaxa){nLoci[5,i]=sum(blah[,i]>125)}
lines(0:(nTaxa-1),nLoci[1,],col="purple")
points(0:(nTaxa-1),nLoci[1,],col="purple",cex=0.5,pch=19)
text(-5,nLoci[1,1]-5,">1500bp",col="purple",cex=0.5)
lines(0:(nTaxa-1),nLoci[2,],col="red")
points(0:(nTaxa-1),nLoci[2,],col="red",cex=0.5,pch=19)
text(-5,nLoci[2,1]-0,">1000bp",col="red",cex=0.5)
lines(0:(nTaxa-1),nLoci[3,],col="orange")
points(0:(nTaxa-1),nLoci[3,],col="orange",cex=0.5,pch=19)
text(-5,nLoci[3,1]+5,">500bp",col="orange",cex=0.5)
lines(0:(nTaxa-1),nLoci[4,],col="green")
points(0:(nTaxa-1),nLoci[4,],col="green",cex=0.5,pch=19)
text(-5,nLoci[4,1]+10,">250bp",col="green",cex=0.5)
lines(0:(nTaxa-1),nLoci[5,],col="cyan")
points(0:(nTaxa-1),nLoci[5,],col="cyan",cex=0.5,pch=19)
text(-5,nLoci[5,1]+15,">125bp",col="cyan",cex=0.5)
arrows(MissingAllowed-1,NLOCI,MissingAllowed-1,NLOCI-10,angle=30,length=0.05)
dev.off()

miss<<-as.matrix(read.table(paste(c("../Results/",project,"_propMissing",folderTail,".txt"),sep="",collapse="")))
miss[miss=="NaN"]=0
pdf(paste(c("../Results/",project,"_AlignmentSummary_PercMissing",folderTail,".pdf"),sep="",collapse=""))
plot(-100,-100,xlim=c(0,(nTaxa-1)),ylim=c(0,max(miss*100)),xlab="Number of Missing Characters Allowed Per Site", ylab="Percent Missing Characters",main=project)
for(i in 1:nTaxa){points(rep(i-1,NLOCI)+runif(NLOCI,-0.25,0.25),miss[,i]*100,cex=0.1)}
points(rep(MissingAllowed-1,NLOCI)+runif(NLOCI,-0.25,0.25),miss[,MissingAllowed]*100,cex=0.1,col="blue")
#arrows(MissingAllowed-1,max(miss*100),MissingAllowed-1,max(miss*100-1),angle=30,length=0.05)
dev.off()

missWeighted=miss*blah
avgMissing=100*apply(missWeighted,2,sum)/apply(blah,2,sum)
pdf(paste(c("../Results/",project,"_AlignmentSummary_AvgPercMissing",folderTail,".pdf"),sep="",collapse=""))
plot(0:(nTaxa-1),avgMissing,xlab="Number of Missing Characters Allowed Per Site",ylab="% of Characters Missing",main=project,pch=19)
lines(0:(nTaxa-1),avgMissing,lwd=2)
points(MissingAllowed-1,avgMissing[MissingAllowed],pch=19,col="blue")
text(MissingAllowed+5,avgMissing[MissingAllowed],paste(floor(avgMissing[MissingAllowed]*100)/100,"%",sep=""),col="blue",cex=0.75)
dev.off()




