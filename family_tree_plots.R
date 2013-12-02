source("family_tree_plotting_functions.R")


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


###plot out family tree with coloring indicating how much they contribute
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



###Picture of your genome across two generations
for(i in 1:10){
	png(file=paste("plots/mother_maternal_grandpars",i,".png",sep=""))
	chr.width=0.3
	chr.lengths<-sapply(1:22,function(chr){max(recoms$mid[recoms$chr==chr])})
	offset=0.2
	par(mar=c(0,0,0,0))
	layout(t(1:3))
	plot.all.chr()
	chr.chunks(my.families[,i],my.col=adjustcolor("purple",0.5),meiosis=1,relly.pos=1); 
	text(0.7,20,"Your genome in\n your Mum",cex=1.5,col="purple")
	plot.all.chr()
	chr.chunks(my.families[,i],my.col=adjustcolor("red",.5),meiosis=2,relly.pos=1); 
	text(0.7,20,paste("Your genome in\n your Maternal\n grandmum"),cex=1.5,col="red")
	plot.all.chr()
	chr.chunks(my.families[,i],my.col=adjustcolor("blue",0.5),meiosis=2,relly.pos=2); 
	text(0.7,20,paste("Your genome \n in your Maternal\n granddad"),cex=1.5,col="blue")
	dev.off()
}

#meiosis<-1; old.relly<-"mother";other.relly<-"sibling's"
#meiosis<-2; old.relly<-"Grandmother";other.relly<-"1st cousin's"
#meiosis<-3; old.relly<-"Great\n Grandmother";other.relly<-"2nd cousin's"
#meiosis<-4; old.relly<-"Great,\n Great Grandmother";other.relly<-"3rd cousin's"
#meiosis<-5; old.relly<-"Great,\n Great, Great\n Grandmother";other.relly<-"4th cousin's"


chr.width=0.3
chr.lengths<-sapply(1:22,function(chr){max(recoms$mid[recoms$chr==chr])})
offset=0.2
par(mar=c(0,0,0,0))
layout(t(1:3))
plot.all.chr()
chr.chunks(my.families[,1],my.col="red",meiosis=meiosis,relly.pos=1); 
text(0.7,20,paste("Your genome in\n your",old.relly),cex=1.5,col="red")
plot.all.chr()
chr.chunks(my.families[,2],my.col=adjustcolor("blue",.5),meiosis=meiosis,relly.pos=1); 
text(0.7,20,paste("Your",other.relly, "\n genome in\n your",old.relly),cex=1.5,col="blue")
plot.all.chr()
chr.chunks(my.families[,1],my.col="red",meiosis=meiosis,relly.pos=1); 
text(0.7,20,paste("Both your genomes \n in your \n",old.relly),cex=1.5,col="purple")
chr.chunks(my.families[,2],my.col=adjustcolor("blue",.5),meiosis=meiosis,relly.pos=1)


#my.families<-replicate(50,simulate.pedigree(num.meioses=10))
#save(file="simulated_pedigrees.Robj",my.families)
load("simulated_pedigrees.Robj")

if(FALSE){
half.rellys<-list()
full.rellys<-list()

for(meiosis in 1:8){
	
	temp<-lapply(1:49,function(i){
		lapply((i+1):50,function(j){
			if(meiosis<6){ 
				IBD.overlap.overall(my.blocks.1=my.families[,i],my.blocks.2=my.families[,j],meiosis=meiosis)
			}else{
				IBD.overlap.overall(my.blocks.1=my.families[,i],my.blocks.2=my.families[,j],meiosis=meiosis,sample.down=TRUE)		
			}

		})
	})
	
	half.rellys[[meiosis]]<-unlist(lapply(temp,function(x){lapply(x,function(y){y$half})}))
	full.rellys[[meiosis]]<-unlist(lapply(temp,function(x){lapply(x,function(y){y$full})}))
	print(meiosis)
	save(file="IBD_sharing_between_rellies.Robj",half.rellys,full.rellys)
}
}
load(file="IBD_sharing_between_rellies.Robj")

png(file="plots/overlap_between_full_cousins.png")
cousin<-c("sibs","first","second","third","fouth","fifth","sixth")
par(mar=c(3,3,2,1))
breaks=seq(-0.5,70,by=1)
layout(matrix(1:6,nrow=2,byrow=TRUE))
for(meiosis in 2:7){
	a<-hist(full.rellys[[meiosis]],breaks=breaks,plot=FALSE)
	upper.x<-a$mid[min(which(a$density==0)) ]
	if(upper.x==0) upper.x<-70
	plot(a$mid,a$density,type="l",xlim=c(0,upper.x))
	points(0:70,dpois(0:70,2*(33.8*(meiosis*2)+22)/(2^(2*meiosis-1))),col="lightgrey",pch=19)
	lines(a$mid,a$density,type="l",lwd=2)
	mtext("probability",2,line=2)
	mtext("# of IBD blocks",1,line=2)
	mtext(paste("full",cousin[meiosis],"cousins"),3,line=0)
}  
dev.off()

png(file="plots/overlap_between_half_cousins.png")
cousin<-c("sibs","first","second","third","fouth","fifth","sixth")
par(mar=c(3,3,2,1))
breaks=seq(-0.5,70,by=1)
layout(matrix(1:6,nrow=2,byrow=TRUE))
for(meiosis in 2:7){
	a<-hist(half.rellys[[meiosis]],breaks=breaks,plot=FALSE)
	upper.x<-a$mid[min(which(a$density==0)) ]
	if(upper.x==0) upper.x<-70
	plot(a$mid,a$density,type="l",xlim=c(0,upper.x))
	points(0:70,dpois(0:70,(33.8*(2*meiosis)+22)/(2^(2*meiosis-1))),col="lightgrey",pch=19)
	lines(a$mid,a$density,type="l",lwd=2)
	mtext("probability",2,line=2)
	mtext("# of IBD blocks",1,line=2)
	mtext(paste("half",cousin[meiosis],"cousins"),3,line=0)
} 
dev.off()

plot(1:8,sapply(half.rellys,function(x){mean(x==0)}),xlab= "ancestor shared k generations back",ylab="Probability of sharing zero genomic regions from anc.",pch=19,cex.lab=1.2)
points(1:8,sapply(full.rellys,function(x){mean(x==0)}),col="red",pch=19)
legend("topleft",legend=c("full-degree relatives","half-degree relatives","approx. full-degree relatives","approx. half-degree relatives"),pch=c(rep(19,2),rep(NA,2)),lty=c(rep(NA,2),rep(1,2)),col=c("red","black","red","black"))
meiosis<-1:8
lines(meiosis,exp(-2*(33.8*(2*meiosis)+22)/(2^(2*meiosis-1))))
lines(meiosis,exp(-(33.8*(2*meiosis)+22)/(2^(2*meiosis-1))),col="red")
