
####generates a transmission chain back through num.meoises following a single lineage
generate.meiotic.path<-function(num.meioses=15,only.1.sex=FALSE){	


	if(only.1.sex %in% c("M","F")){ 	
		my.path<-sample(inds.sex[grep(only.1.sex,inds.sex)],num.meioses)
		}else{
		my.path<-sample(inds.sex,num.meioses)
	}
	
	chunks<-sapply(1:22,function(chr){
		chunks<-vector("list",num.meioses)

		chr.length<-max(recoms$mid[recoms$chr==chr])
		my.chunks<-data.frame(l=0,r=chr.length+1)
		
	
		for(meiosis in 1:num.meioses){
			chunks[[meiosis]]<-my.chunks

		#	print(my.chunks)
			new.chunks<-break.chrom(ind.sex=my.path[meiosis],chr=chr,my.chunks)
			if(length(new.chunks)==0){break;}
			my.chunks<-new.chunks
			stopifnot(ncol(new.chunks)==2 | ncol(new.chunks)==0)
		#	cat("my new chunks\n")
		#	print(new.chunks)
		#recover()
		}	
	chunks
	})
}


##initiate all the chr for simulate.pedigree
initiate.chrs<-function(){
my.chunks.all.chr<-sapply(1:22,function(chr){
	chr.length<-max(recoms$mid[recoms$chr==chr])
	my.chunks<-data.frame(l=0,r=chr.length+1)
	return(list(my.chunks))
})
my.chunks.all.chr
}


###simulate transmission through the entire pedigree
simulate.pedigree<-function(num.meioses){
 
	family.chunks<-sapply(1:(num.meioses+1),function(meiosis){vector("list",2^(meiosis))})
	family.chunks[[1]][[1]]<-initiate.chrs()
	family.chunks[[1]][[2]]<-initiate.chrs()
	
	for(meiosis in 1:num.meioses){
		par.sex<-rep(c("M","F"),2^(meiosis-1))   ##mother father
		
		for(my.par in 1:2^meiosis){
	
			ind.sex<-sample(inds.sex[grep(par.sex[my.par],inds.sex)],1)
			my.chunks.ind<-family.chunks[[meiosis]][[my.par]]
			chunks<-lapply(1:22,function(chr){
					my.chunks<-my.chunks.ind[[chr]]
					if(length(my.chunks)==0) return(NULL)
					new.chunks<-break.chrom.2(ind.sex=ind.sex,chr=chr,my.chunks)
		#			recover()
		#			if(length(new.chunks)==0){return("blah");}
			#		stopifnot(ncol(new.chunks)==2 | ncol(new.chunks)==0)
					new.chunks
				})	
			g.par1<-lapply(chunks,function(x){x[[1]]})
			g.par2<-lapply(chunks,function(x){x[[2]]})
			family.chunks[[meiosis+1]][[(my.par-1)*2+1]] <- g.par1
			family.chunks[[meiosis+1]][[(my.par-1)*2+2]] <- g.par2
			cat((my.par-1)*2+1,(my.par-1)*2+2," ")
		}
	}
family.chunks
}


####simulates transmission back through meiosis for a chromosome breaking the chunks that it is passed at recombination events in the individual supplied
break.chrom<-function(ind.sex,chr,my.chunks){
	stopifnot(ncol(my.chunks)==2 | ncol(my.chunks)==0)
	pos<-recoms$mid[recoms$chr==chr & recoms$sexind==ind.sex]
	pos<-sort(pos)
	if(length(pos)==0){ ##no rec in chunk this gen.
			 if(sample(0:1,1)){ 
			 	return(data.frame(my.chunks))
			 } else{
			 	return(NULL)
			 }
	}		
	new.chunks<-apply(my.chunks,1,function(chunk){
		my.recs<- pos[chunk["l"] < pos & pos < chunk["r"]]
	
		if(length(my.recs)==0){ ##no rec in chunk this gen.
			 if(sample(0:1,1)){ 
			 	return(data.frame(l=chunk["l"], r=chunk["r"]))
			 	} else{
			 		return(NULL)
			 	}
			 }
		intervals<- c(chunk["l"] ,my.recs, chunk["r"]) #break up chunk by rec
		odds<-seq(1,length(intervals)-1,by=2)
		evens<-seq(2,length(intervals)-1,by=2)
		par1<-data.frame(cbind(l=intervals[odds],r=intervals[odds+1]))   ## left and right boundaries of odd segments 
		par2<-data.frame(cbind(l=intervals[evens],r=intervals[evens+1]))  ## left and right boundaries of even segments 
#		recover()
		if(sample(0:1,1)){ 
			return(par1)
			}else{
			return(par2)	
				}
				##commented out text to return both sets of segments
#		if(sample(0,1)){ 
#			return(data.frame(cbind(par1=par1,par2=par2)))
#		}else{
#			return(data.frame(cbind(par1=par2,par2=par1)))	
#		}
	})	
#	recover()

	my.new.chunks<-data.frame(l=c(unlist(sapply(new.chunks,function(x){x$l}))),r=	c(unlist(sapply(new.chunks,function(x){x$r}))))	

	stopifnot(all(my.new.chunks$l <= my.new.chunks$r))
	stopifnot(ncol(my.new.chunks)==2 | ncol(my.new.chunks)==0)
	
#		cat("my new chunks in func\n")
#	print(my.new.chunks)

	return(my.new.chunks)
	}


####simulates transmission back through meiosis for a chromosome breaking the chunks that it is passed at recombination events in the individual supplied. This function keeps track of both grandparental contributions so that entire pedigrees can be simulated

break.chrom.2<-function(ind.sex,chr,my.chunks){
	stopifnot(ncol(my.chunks)==2 | ncol(my.chunks)==0)
	pos<-recoms$mid[recoms$chr==chr & recoms$sexind==ind.sex]
	pos<-sort(pos)
	if(length(pos)==0){ ##no rec in chunk this gen.
			 if(sample(0:1,1)){ 
			 	return(list(data.frame(my.chunks),NULL))
			 } else{
			 	return(list(NULL,data.frame(my.chunks)))
			 }
	}		
	new.chunks<-apply(my.chunks,1,function(chunk){
		my.recs<- pos[chunk["l"] < pos & pos < chunk["r"]]
	
		if(length(my.recs)==0){ ##no rec in chunk this gen.
			 if(sample(0:1,1)){ 
			 		return(list(data.frame(l=chunk["l"], r=chunk["r"]),NULL))
			 	} else{
			 		return(list(NULL,data.frame(l=chunk["l"], r=chunk["r"])))
			 	}
			 }
		intervals<- c(chunk["l"] ,my.recs, chunk["r"]) #break up chunk by rec
		odds<-seq(1,length(intervals)-1,by=2)
		evens<-seq(2,length(intervals)-1,by=2)
		par1<-data.frame(cbind(l=intervals[odds],r=intervals[odds+1]))   ## left and right boundaries of odd segments 
		par2<-data.frame(cbind(l=intervals[evens],r=intervals[evens+1]))  ## left and right boundaries of even segments 

#				recover()
		if(sample(0:1,1)){ 
			return(list(par1,par2))
		}else{
			return(list(par2,par1))	
		}
	})	
#	recover()

	par1<-lapply(new.chunks,function(x){x[[1]]})
	par2<-lapply(new.chunks,function(x){x[[2]]})
	chunks.par1<-data.frame(l=c(unlist(sapply(par1,function(x){x$l}))),r=	c(unlist(sapply(par1,function(x){x$r}))))	
	chunks.par2<-data.frame(l=c(unlist(sapply(par2,function(x){x$l}))),r=	c(unlist(sapply(par2,function(x){x$r}))))	
	stopifnot(all(chunks.par1$l <= chunks.par1$r))
	stopifnot(all(chunks.par2$l <= chunks.par2$r))
stopifnot((ncol(chunks.par1)==2 | ncol(chunks.par1)==0) & (ncol(chunks.par2)==2 | ncol(chunks.par2)==0))

	if(ncol(chunks.par1)==0) chunks.par1<-NULL
	if(ncol(chunks.par2)==0) chunks.par2<-NULL
	
		
#		cat("my new chunks in func\n")
#	print(my.new.chunks)

	return(list(chunks.par1,chunks.par2))
	}


###Functions to process blocks
	
get.num.blocks<-function(all.chunks){
	num.blocks.all<-sapply(1:15,function(meiosis){
		num.blocks<-sapply(1:dim(all.chunks)[3],function(sim){
		sum(unlist(lapply(all.chunks[meiosis,,sim],function(x){if(is.null(x)){return(0)}; nrow(x)})))
		})
	})
num.blocks.all
}

get.length.blocks<-function(all.chunks){
	blocks.length.all<-sapply(1:15,function(meiosis){
		num.blocks<-sapply(1:dim(all.chunks)[3],function(sim){
		sum(unlist(lapply(all.chunks[meiosis,,sim],function(x){if(is.null(x)){return(0)}; sum(x[,2]-x[,1])})))
		})
	})
blocks.length.all
}

num.zeros<-function(num.blocks.all){
sapply(1:15,function(meiosis){
	mean(num.blocks.all[,meiosis]==0)
	})
}


get.prop.zero.chr<-function(all.chunks,chr){
	num.blocks.chr<-sapply(1:15,function(meiosis){
		num.blocks<-sapply(1:dim(all.chunks)[3],function(sim){
		lapply(all.chunks[meiosis,chr,sim],function(x){if(is.null(x)){return(0)}; nrow(x)})
		})
	})

apply(num.blocks.chr,2,function(x){mean(x==0)})
}
