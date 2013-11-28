
####Simple family tree
png(file="family_tree.png",width = 800, height = 800)
k=8
pie<-lapply(1:2^k,function(x){c(k,k+1)})
repeat.cols<-c("red","blue") #c("white","white") #
my.cols<-rep(repeat.cols,2^(k-1))


my.pie(pie,sector.colors=my.cols)

for(i in (k-1):1){
	pie<-lapply(1:2^i,function(x){c(i,i+1)})
	my.cols<-rep(repeat.cols,2^(i-1))

	my.pie(pie,add=TRUE,sector.colors=my.cols)
}
dev.off()

##simulate a fmaily running back num.meiosis generations
family.chunks<-simulate.pedigree(num.meioses=11)

pdf(file="family_tree_w_11_gens.pdf",width = 800, height = 800)
k=11
par(mar=c(0,0,0,0))
pie<-lapply(1:2^k,function(x){c(k,k+1)})
repeat.cols<-c("red","blue") #c("white","white") #
my.cols<-rep(repeat.cols,2^(k-1))

num.blocks<-sapply(family.chunks[[k]],function(ind){ sum(unlist(lapply(ind,function(x){if(is.null(x)){return(0)}; nrow(x)})))})
amount.genome<-unlist(lapply(family.chunks[[k]],function(ind){ sum(unlist(lapply(ind,function(x){if(is.null(x)){return(0)}; sum(x[,2]-x[,1])})))}))
frac.genome<-amount.genome/genome.length

#my.cols<-sapply(1:length(my.cols),function(i){adjustcolor(my.cols[i],0.03+0.97*frac.genome[i])})
my.cols<-sapply(1:length(my.cols),function(i){adjustcolor(my.cols[i],frac.genome[i]/max(frac.genome))})
my.cols[num.blocks==0]<-"white"

my.pie(pie,sector.colors=my.cols)
	text(-k-0.5,-0.6,paste(format(100*mean(frac.genome)/2,digit=2), format(100*mean(frac.genome==0),digit=2),sep=", "),cex=1,srt=90)
	text(k+0.5,-0.6,paste(format(min(100*frac.genome)/2,digit=2),"-",format(100*max(frac.genome)/2,digit=2)),cex=1,srt=90)
#	text(-k/2-1,-1.5,"Mean % contribution, % with zero",cex=1)
#text(k/2+1,-1.5,"Min-Max % contribution",cex=1)

for(i in (k-1):1){
	pie<-lapply(1:2^i,function(x){c(i,i+1)})
	my.cols<-rep(repeat.cols,2^(i-1))
	num.blocks<-sapply(family.chunks[[i]],function(ind){ sum(unlist(lapply(ind,function(x){if(is.null(x)){return(0)}; nrow(x)})))})
	amount.genome<-unlist(lapply(family.chunks[[i]],function(ind){ sum(unlist(lapply(ind,function(x){if(is.null(x)){return(0)}; sum(x[,2]-x[,1])})))}))
	frac.genome<-amount.genome/genome.length
	text(-i-0.5,-0.6,paste(format(100*mean(frac.genome)/2,digit=2), format(100*mean(frac.genome==0),digit=2),sep=", "),cex=1,srt=90)
	text(i+0.5,-0.6,paste(format(100*min(frac.genome)/2,digit=2),"-",format(100*max(frac.genome)/2,digit=2)),cex=1,srt=90)
#	my.cols<-sapply(1:length(my.cols),function(i){adjustcolor(my.cols[i],0.03+0.97*frac.genome[i])})
my.cols<-sapply(1:length(my.cols),function(i){adjustcolor(my.cols[i],frac.genome[i]/max(frac.genome))})
	my.cols[num.blocks==0]<-"white"

	my.pie(pie,add=TRUE,sector.colors=my.cols)
}

dev.off()



######chr blocks for individuals
png(file="family_tree_w_trans_of_chr_chunks_4_gens.png")
k=4; chr=1
pie<-lapply(1:2^k,function(x){c(k,k+1)})
num.anc<-2^k; 

num.blocks<-sapply(family.chunks[[k]],function(ind){ if(is.null(ind[[chr]])){return(0)};   nrow(ind[[chr]])})

anc.cols<-rep(NA,num.anc)
 #anc.cols<-rainbow(num.anc) 
anc.cols[num.blocks==0]<-"white"
anc.cols[num.blocks!=0]<- (rainbow(sum(num.blocks!=0)))

my.pie(pie,sector.colors=anc.cols,lower.y=-6)
	these.ancs.mum<-1:2^(k-1)
	these.ancs.dad<-(2^(k-1)+1):2^k

sapply(these.ancs.mum,plot.blocks,meiosis=k,which.side="left")
sapply(these.ancs.dad,plot.blocks,meiosis=k,which.side="right")
text(-k/2,y=-0.5,"chromosome from mum\n over generations")
text(k/2,y=-0.5,"chromosome from dad\n over generations")
for(i in (k-1):1){
	pie<-lapply(1:2^i,function(x){c(i,i+1)})
	num.anc<-2^i; 
	num.blocks<-sapply(family.chunks[[i]],function(ind){ if(is.null(ind[[chr]])){return(0)};   nrow(ind[[chr]])})
	anc.cols<-rep(NA,num.anc)
	 #anc.cols<-rainbow(num.anc) 
	anc.cols[num.blocks==0]<-"white"
	anc.cols[num.blocks!=0]<- (rainbow(sum(num.blocks!=0)))

	my.pie(pie,add=TRUE,sector.colors=anc.cols)
	these.ancs.mum<-1:2^(i-1)
	these.ancs.dad<-(2^(i-1)+1):2^i

sapply(these.ancs.mum,plot.blocks,meiosis=i,which.side="left")
sapply(these.ancs.dad,plot.blocks,meiosis=i,which.side="right")

}

dev.off()



#which.relly<-2; old.relly<-"mother";other.relly<-"sibling's"
#which.relly<-3; old.relly<-"Grandmother";other.relly<-"1st cousin's"
#which.relly<-4; old.relly<-"Great\n Grandmother";other.relly<-"2nd cousin's"
#which.relly<-5; old.relly<-"Great,\n Great Grandmother";other.relly<-"3rd cousin's"
which.relly<-6; old.relly<-"Great,\n Great, Great\n Grandmother";other.relly<-"4th cousin's"


chr.width=0.3
chr.lengths<-sapply(1:22,function(chr){max(recoms$mid[recoms$chr==chr])})
offset=0.2
par(mar=c(0,0,0,0))
layout(t(1:3))
plot.all.chr()
chr.chunks(family.1.chunks[[which.relly]],my.col="red"); 
text(0.7,20,paste("Your genome in\n your",old.relly),cex=1.5,col="red")
plot.all.chr()
chr.chunks(family.2.chunks[[which.relly]],my.col=adjustcolor("blue",.5)); 
text(0.7,20,paste("Your",other.relly, "\n genome in\n your",old.relly),cex=1.5,col="blue")
plot.all.chr()
chr.chunks(family.1.chunks[[which.relly]],my.col="red"); 
text(0.7,20,paste("Both your genomes \n in your \n",old.relly),cex=1.5,col="purple")
chr.chunks(family.2.chunks[[which.relly]],my.col=adjustcolor("blue",.5))





sapply(1:22,function(chr){  
   	blocks<-family.2.chunks[[3]][[1]][[chr]]/chr.lengths[1] 
	if(!is.null(nrow(blocks))) apply(blocks,1,function(x){x= polygon(x=c(x,rev(x)),c(rep(chr+offset,2),rep(chr+offset+chr.width,2)),col=adjustcolor("blue",.5))})
	 blocks<-family.2.chunks[[3]][[2]][[chr]]/chr.lengths[1]
 	if(!is.null(nrow(blocks))) apply(blocks,1,function(x){x= polygon(x=c(x,rev(x)),c(rep(chr-offset,2),rep(chr-offset+chr.width,2)),col=adjustcolor("blue",.5))})
})

  
