
recoms<-read.table("recombination_events.out",as.is=TRUE,head=TRUE)
source("transmission_sims_functions.R")

inds.sex<-unique(recoms$sexind)

#####Probability of transmitting all chromosomes without recombination, for mothers "M", change to "F" for fathers.
inds<-unique(temp1$ind)
prob.no.rec<-sapply(1:22,function(chr){
	my.chr<-recoms[recoms$par=="M" & recoms$chr==chr,]
	mean(!(inds %in% my.chr$ind))
	})
	
	prod(prob.no.rec)
	
	
	chr.transmitted<-function(ind,chr,sex,chr.length){
	pos<-recoms$mid[recoms$chr==chr & recoms$par==sex & recoms$ind==ind]
	intervals<-diff(c(0,sort(pos),chr.length))
	if(length(intervals)==1) return(c(0,intervals))
	odds<-seq(1,length(intervals),by=2)
	evens<-seq(2,length(intervals),by=2)
	return(c(sum(intervals[odds]),sum(intervals[evens])))
	}
	
	
###find out how much grandparental material both fathers and mothers transmit.	
trans.male<-list()
trans.female<-list()
for(chr in 1:22){
	chr.length<-max(recoms$mid[recoms$chr==chr])
	chr.trans.male<-sapply(inds,chr.transmitted,sex="F",chr=chr,chr.length=chr.length)
	chr.trans.female<-sapply(inds,chr.transmitted,sex="M",chr=chr,chr.length=chr.length)

	trans.male[[chr]]<-chr.trans.male
	trans.female[[chr]]<-chr.trans.female
	}

###transmission of grandparental material both fathers and mothers per chromosome
layout(1:2)
pdf(file="dist_transmitted_chromosome.pdf")	
	for(chr in 1:22){
	chr.length<-max(recoms$mid[recoms$chr==chr])
	layout(c(1:2))
	chr.trans.male<-trans.male[[chr]]
	chr.trans.female<-trans.female[[chr]]
	
	hist(c(chr.trans.male),main=paste("male chr. ",chr),xlab="amount transmitted",breaks=seq(0,chr.length,length=20))
	hist(c(chr.trans.female),main=paste("female chr. ",chr),xlab="amount transmitted",breaks=seq(0,chr.length,length=20))
}
dev.off()


##reformat list of transmitted material per chr. per ind to be an array
tot.trans.male<-array(dim=c(length(inds),22,2))
tot.trans.female<-array(dim=c(length(inds),22,2))

for(chr in 1:22){
	tot.trans.female[,chr,]<-trans.female[[chr]]
	tot.trans.male[,chr,]<-trans.male[[chr]]

	}


###transmission of grandparental material both fathers and mothers simulated for entire genome
genome.length<-sum(sapply(1:22,function(chr){max(recoms$mid[recoms$chr==chr])}))
pdf("fraction_of_genome_transmitted.pdf")
layout(1:2)
tot.genome.trans<-sapply(1:length(inds),function(ind){replicate(100,sum(sapply(1:22,function(chr){tot.trans.female[ind,chr,sample(2,1)]})))})   ##for each individual simulate whether they transmit grandmat. or grandpat material on each chr
hist(tot.genome.trans/genome.length,main="fraction of grandparental genome transmitted via females",xlab="fraction")
tot.genome.trans<-sapply(1:length(inds),function(ind){replicate(100,sum(sapply(1:22,function(chr){tot.trans.male[ind,chr,sample(2,1)]})))})
hist(tot.genome.trans/genome.length,main="fraction of grandparental genome transmitted via males",xlab="fraction")
dev.off()



###generate 5000 meiotic paths of transmission 
##this takes a while, so it is commented out.
if(FALSE){
all.chunks<-replicate(5000,generate.meiotic.path)
all.chunks.fem<-replicate(5000,generate.meiotic.path(only.1.sex="M"))
all.chunks.mal<-replicate(5000,generate.meiotic.path(only.1.sex="F"))
save(file="my.chunks.Robj",all.chunks,all.chunks.fem,all.chunks.mal)
}
##you can load the Rboj instead
load("my.chunks.Robj")

	
##count up all of the blocks inherited from a specifc ancestor
num.blocks.all<-get.num.blocks(all.chunks)
num.blocks.fem<-get.num.blocks(all.chunks.fem)
num.blocks.mal<-get.num.blocks(all.chunks.mal)

###Figures number of blocks inherited from a specifc ancestor
these.meioses<- c(2,4,7,9) # 1:15
png(file="distribution_num_blocks_egs.png")
layout(matrix(1:4,nrow=2,byrow=TRUE))
add.legends<-rep(TRUE,15)
add.legends[these.meioses]<-c(rep(FALSE,3),TRUE)
#pdf(file="~/Dropbox/Pedigree_Recom/distribution_num_blocks.pdf")
sapply(these.meioses,function(meiosis){
	tmp<-hist(num.blocks.mal[,meiosis],breaks=seq(-.9,50+.1,by=1),freq=FALSE,plot=FALSE)
	plot(tmp$mids,tmp$intensities,type="l",col="blue",xlab="Number of blocks",ylab="Probability of having # blocks",main=paste(meiosis,"generations back"),lwd=2)
	if(add.legends[meiosis]) legend("topright",c("typical transmissions","Matrilineal","Patrilineal"),col=c("black","red","blue"),lty=2,lwd=2)
	tmp<-hist(num.blocks.all[,meiosis],breaks=seq(-.9,50+.1,by=1),freq=FALSE,plot=FALSE)
	lines(tmp$mids,tmp$intensities,type="l",lwd=2)
	tmp<-hist(num.blocks.fem[,meiosis],breaks=seq(-.9,50+.1,by=1),freq=FALSE,plot=FALSE)
	lines(tmp$mids,tmp$intensities,type="l",col="red",lwd=2)
		
#	points(0:45-.1,dpois(0:45,(33*(meiosis-1)+22)/(2^(meiosis-1))),col="blue")
	})
dev.off()

##compare to analytical approximation, based on poisson distribution of number of blocks
these.meioses<-  1:15
png(file="~/Dropbox/Pedigree_Recom/distribution_num_blocks_vs_approx.png")
layout(matrix(1:4,nrow=2,byrow=TRUE))
add.legends<-rep(TRUE,15)
#add.legends[these.meioses]<-c(rep(FALSE,3),TRUE)
pdf(file="~/Dropbox/Pedigree_Recom/distribution_num_blocks_vs_approx.pdf")
sapply(these.meioses,function(meiosis){
	tmp<-hist(num.blocks.all[,meiosis],breaks=seq(-.9,50+.1,by=1),freq=FALSE,plot=FALSE)
	plot(tmp$mids,tmp$intensities,type="l",col="blue",xlab="Number of blocks",ylab="Probability of having # blocks",main=paste(meiosis,"generations back"),lwd=2)
	if(add.legends[meiosis]) legend("topright",c("typical transmissions","Matrilineal","Patrilineal"),col=c("black","red","blue"),lty=2,lwd=2)

	points(0:45-.1,dpois(0:45,(33*(meiosis-1)+22)/(2^(meiosis-1))),col="blue")
	})
dev.off()



##sum up total genomic length of the blocks inherited from a specifc ancestor
length.blocks.all<-get.length.blocks(all.chunks)
length.blocks.mal<-get.length.blocks(all.chunks.mal)
length.blocks.fem<-get.length.blocks(all.chunks.fem)

##graph distribution of total genomic length of the blocks inherited from a specifc ancestor
these.meioses<- c(2,4,7,9) # 1:15
png(file="distribution_amount_autosomes_egs.png")
layout(matrix(1:4,nrow=2,byrow=TRUE))
add.legends<-rep(TRUE,15)
add.legends[these.meioses]<-c(rep(FALSE,3),TRUE)
#pdf(file="~/Dropbox/Pedigree_Recom/distribution_amount_autosomes.pdf")
sapply(these.meioses,function(meiosis){
	my.range<-range(c(length.blocks.fem[,meiosis],length.blocks.mal[,meiosis],length.blocks.all[,meiosis]) /genome.length)
	my.breaks<-seq(my.range[1],my.range[2],length=20)
	tmp<-hist(length.blocks.fem[,meiosis]/genome.length,freq=FALSE,plot=FALSE,breaks=my.breaks)
	plot(tmp$mids,tmp$intensities,type="l",col="red",xlab="Fraction of autosome",ylab="Probability",main=paste(meiosis,"generations back"),lwd=2)
	if(add.legends[meiosis]) legend("topright",c("typical transmissions","Matrilineal","Patrilineal"),col=c("black","red","blue"),lty=2,lwd=2)
	tmp<-hist(length.blocks.mal[,meiosis]/genome.length,freq=FALSE,plot=FALSE,breaks=my.breaks)
	lines(tmp$mids,tmp$intensities,lwd=2,col="blue")
	tmp<-hist(length.blocks.all[,meiosis]/genome.length,freq=FALSE,plot=FALSE,breaks=my.breaks)
	lines(tmp$mids,tmp$intensities,col="black",lwd=2)
		
#	points(0:45-.1,dpois(0:45,(33*(meiosis-1)+22)/(2^(meiosis-1))),col="blue")
	})
dev.off()

###Probability of receiving zero of your genome from an ancestor k meioses back
png(file="Prob_zero_blocks.png")
plot(1:15,num.zeros(num.blocks.all),type="b",xlab="Generations back",ylab="Probability of inheriting zero blocks of ancestral genome")
lines(1:15,num.zeros(num.blocks.fem),type="b",col="red")
lines(1:15,num.zeros(num.blocks.mal),type="b",col="blue")
legend("topleft",c("typical transmissions","Matrilineal","Patrilineal"),col=c("black","red","blue"),lty=1,lwd=2)
dev.off()

length.blocks.all<-sapply(1:10,function(meiosis){
	num.blocks<-sapply(1:dim(all.chunks)[3],function(sim){
	unlist(lapply(all.chunks[2,,1],function(x){x[,2]-x[,1]}))
	})
})


##prob. of zero blocks inherited for each chromosome
prop.zero.chr<-sapply(1:22,function(chr){get.prop.zero.chr(all.chunks,chr=chr)})	

png(file="Prob_zero_blocks_chr.png")
plot(c(1,15),c(0,1),type="n",ylab="Prob. of inheriting zero blocks of ancestral chromosome",
xlab="generations back", main="Probability of inheriting zero blocks from a particular ancestor")		
my.cols<-rainbow(22,end =0.75);sapply(1:22,function(i){lines(prop.zero.chr[,i],col=my.cols[i])})
text(10,0.45,"chromosome")
text(seq(8,14,length=11),rep(0.4,11),1:11,col=my.cols[1:11])
text(seq(8,14,length=11),rep(0.35,11),12:22,col=my.cols[12:22])
dev.off()


#####Compare probability of zero from sims to theoretical approximation
k<-1:15
png(file="Prob_zero_blocks_vs_theory.png")
plot(k,1-exp(-(22+35*(k-1))/(2^(k-1))),xlab="generation back",ylab="Probability that an ancestor is an autosomal genetic ancestor",type="l");points(k,1-num.zeros(num.blocks.all))
 legend("bottomleft",legend=c("Simulations","Theoretical approximation"),pch=c(1,NA),lty=c(NA,1))
dev.off()


####Number of genetic ancestors vs genealogical ancestors
png(file="Num_genetics_vs_genealogical_ancs.png")
k<-1:20;num.blocks<-(22+33*(k-1))/(2^(k-1));

plot(2*(1-exp(-num.blocks))*(2^(k-1)),xlab="Generations ago",ylab="Expected number of ancestors");lines(2^k,col="red")
legend(x="topleft",legend=c("Genealogical","Genetic"),lty=c(1,NA),pch=c(NA,1),col=c("red","black"))
dev.off()








	

		
		