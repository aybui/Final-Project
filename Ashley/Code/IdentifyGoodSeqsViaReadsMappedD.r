identifyGoodSeqs<<-function(begInd,endInd,nHomologs,minReadsMapped){

	firstWrite=T
	for(ind in begInd:endInd){

		if(!file.exists(paste("../I",ind,"/I",ind,"_nMapped.txt",collapse="",sep=""))){next}

		for(i in 1:nHomologs){
		
			blah<<-read.table(paste("../I",ind,"/I",ind,"_nMapped.txt",collapse="",sep=""))
			goodLoci=blah[blah[,3]>=minReadsMapped & blah[,2]==i,1]
			
			towrite=cbind(rep(i,length(goodLoci)),rep(ind,length(goodLoci)),goodLoci)
			if(length(towrite)>0){
				write(t(towrite),ncol=3,sep="\t",file="../Results/GoodConSeqIDs.txt",append=!firstWrite)
				firstWrite=F
			}			

			if(i==1){		nMapped1<<-blah[blah[,2]==i,3]
			}else if(i==2){	nMapped2<<-blah[blah[,2]==i,3]
			}else if(i==3){	nMapped3<<-blah[blah[,2]==i,3]
			}else if(i==4){	nMapped4<<-blah[blah[,2]==i,3]
			}else if(i==5){	nMapped5<<-blah[blah[,2]==i,3]
			}else if(i==6){	nMapped6<<-blah[blah[,2]==i,3]
			}
			
		}
		
		write(as.integer(round(quantile(nMapped1,0.05,names=FALSE))),file="minReads.txt",append=TRUE)
	}

}

args <- commandArgs(TRUE)
if(file.exists("minReads.txt")){file.remove("minReads.txt")}
identifyGoodSeqs(args[1], args[2], args[3], args[4])
