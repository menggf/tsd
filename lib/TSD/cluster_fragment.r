clustervector<-function(x, rg, min.split.reads=1){
	cl=rep(0,length(x));
	zz=0;
	center=vector();
	while(length(cl[cl==0]) > 0){
		zz=zz+1;
		tx=table(x[cl==0]);
		cc=as.numeric(names(which.max(tx)))
		center=append(center, cc);
		tag=cl==0 & x > cc-rg & x < cc+rg;
    if(length(tag[tag]) < min.split.reads) 
      break;
		cl[tag]=zz;
	}
	names(center)=1:zz
	return(list(cl=cl, center=center))
}
args <- commandArgs(TRUE)

da<-read.table(paste(args[1],"/sites.txt",sep=""), sep="\t", fill=TRUE)
da=da[!is.na(as.vector(da$V1)) & as.vector(da$V1) != "" & !is.na(as.vector(da$V2)) & !is.na(as.vector(da$V3))  ,]
chr=as.vector(da$V1)
from=as.numeric(as.vector(da$V2));
to=as.numeric(as.vector(da$V3));
tag=!is.na(from) & !is.na(to)
chr=chr[tag]
to=to[tag]
from=from[tag]


cc=0;
labs=rep(0, length(chr))
poss1=rep(-1, length(chr))
poss2=rep(-1, length(chr))
for(tg in unique(chr)){
   #print(tg)
  wh=chr==tg;
  cl1=vector()
  cl2=vector()
  if(tg != "insert"){
    cl1=clustervector(from[wh], 200, 1);
    cl2=clustervector(to[wh],   200, 1);
  }else{
    cl1=clustervector(from[wh], 20, 1);
    cl2=clustervector(to[wh],   20, 1);
  }
  labs[wh]=paste(cl1$cl+cc, cl2$cl +cc, sep="#");
  cc=cc+length(cl1$cl);
  poss1[wh]=cl1$center[cl1$cl]
  poss2[wh]=cl2$center[cl2$cl]
}

new.da=cbind(da[tag,], site=paste("sites",labs,sep=""), pos1=poss1, pos2=poss2)
write.table(new.da, paste(args[1],"/new_sites.txt",sep=""),row.names=F,col.names=F,sep="\t",quote=F)


