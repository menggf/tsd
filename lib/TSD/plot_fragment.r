args <- commandArgs(TRUE)
#print(args)
#cv=as.vector(args[7:length(args)])
#left=args[4]
#ns=args[5]
#right=args[6]
#print(args)

input<-scan(args[2],sep="\n", what=list(""), quiet=TRUE)[[1]]
nn=length(input)
s0=strsplit(input[1],"\t")[[1]];
start = as.numeric(s0[3]);
end = as.numeric(s0[4]);
s=strsplit(input[2],"\t")[[1]];

nfrag=(length(s)+2)/6;

#used=seq_len(nfrag)
#n.used=length(used)
seps=as.numeric(strsplit(args[5],",")[[1]])
seps2=strsplit(args[6],",")[[1]]
pdf(paste(args[1],"/seq_pic/seq",args[3],".pdf",sep=""), height=max(4, (nn-3)/1.5),width=max(7,1.5*nfrag))
plot(c((start)*6-1,(end+1)*6+2),c(1,-2*(nn-2)), type="n", xlab="",ylab="",xaxt='n', yaxt='n', main=args[4])
mycols=vector()
for(i in start:end+1){
	if(s[(i-1)*6+2] =="1"){
		arrows((i-1)*6+0.5, 0, (i-1)*6+5.5,  0, angle = 5, lwd=4, col="blue")
		text(x=(i-1)*6+2.5, y=0.2, labels=paste(s[(i-1)*6+1],":", s[(i-1)*6+3],"-",s[(i-1)*6+4],sep=""),pos=3, cex=0.6)
		mycols[i]=colors()[89];
	}else{
		arrows((i-1)*6+5.5, 0, (i-1)*6+0.5,  0 , angle = 5, lwd=4, col="green")
		text(x=(i-1)*6+2.5, y=0.2, labels=paste(s[(i-1)*6+1],":", s[(i-1)*6+3],"-",s[(i-1)*6+4],sep=""),pos=3, cex=0.6)
		mycols[i]=colors()[90];
	}
}
if(!is.na(seps[1])){
	for(i in seq_along(seps)){
		if(seps2[i] == "l"){
			points(x=seps[i]*6, y=0, pch="<", col="red")
		}else{
			points(x=seps[i]*6, y=0, pch=">", col="red")
		}
	}
}
for(i in 3:nn){
	ss=as.numeric(strsplit(input[i],"\t")[[1]]);
	mm=length(ss);
	ms=(mm-1)/2;
	for(j in 1:((mm-1)/2)+1){
		if(ss[j+ms]==-1)
			next();
		lines(c(ss[j]*6+0.5,ss[j]*6+5.5), c(-2 *(i-3)-1.2, -2 *(i-3)-1.2), lwd= 4, col=mycols[ss[j]+1])
		
		if(ss[2+ms]!=0 & ss[2+ms]!=-1& ss[2]==0){
			text(0, -2 *(i-3)-1.2, labels="<",col="blue")
		}
		#if(ss[mm] == nfrag-1 and ss[ms+ms+1]!=){
		#	text(nfrag*6, -2 *(i-3)-1.2, labels=">",col="green")
		#}
		text(end*6+7.5, -2 *(i-3)-1.2, labels=ss[1],col="brown")
	}
}
dev.off()
