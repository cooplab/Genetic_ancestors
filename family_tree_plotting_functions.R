plot.blocks<-function(this.anc,meiosis,chr=1,which.side){	
	blocks<-family.chunks[[meiosis]][[this.anc]][[chr]];
	if(length(blocks)==0){ return()}
	chr.length<-max(recoms$mid[recoms$chr==chr])
	if(which.side=="left") anc.pos<- -1*meiosis-0.5
	if(which.side=="right")  anc.pos<- meiosis + 0.5 
blocks<-(blocks/chr.length)*(-5)-1
#	segments(anc.pos,blocks$l,anc.pos,blocks$r,col=anc.cols[this.anc],lwd=3,lend=2)
	apply(blocks,1,function(x){x= polygon(c(rep(anc.pos,2),rep(anc.pos+.4,2)) , y=c(x,rev(x)),col=anc.cols[this.anc],border="lightgrey")})
}


##draw paternal and maternal chrs coloured by where IBD comes from, meiosis generations back, in ancestor relly.pos
chr.chunks<-function(my.blocks,my.col,meiosis,relly.pos){
	num.ancs<-1:2^(meiosis+1)
	maternal.chr<-num.ancs[(num.ancs %% 2)==0][relly.pos]
	paternal.chr<-num.ancs[(num.ancs %% 2)==1][relly.pos]

	sapply(1:22,function(chr){  
	   	blocks<-my.blocks[[meiosis+1]][[maternal.chr]][[chr]]/chr.lengths[1] 
		if(!is.null(nrow(blocks))) apply(blocks,1,function(x){x= polygon(x=c(x,rev(x)),c(rep(chr+offset,2),rep(chr+offset+chr.width,2)),col=my.col)})
		blocks<-my.blocks[[meiosis+1]][[paternal.chr]][[chr]]/chr.lengths[1]
	 	if(!is.null(nrow(blocks))) apply(blocks,1,function(x){x= polygon(x=c(x,rev(x)),c(rep(chr-offset,2),rep(chr-offset+chr.width,2)),col=my.col)})
	})
}


plot.all.chr<-function(){
plot(c(0,1+.03),c(1,22),type="n",axes=FALSE,xlab="",ylab="")
box()
sapply(1:22,function(chr){
	scale.chr<-c(0,chr.lengths[chr]/chr.lengths[1])
	text(chr.lengths[chr]/chr.lengths[1]+.05,chr+.25,chr)
	polygon(x=c(scale.chr,rev(scale.chr)),y=c(rep(chr+offset,2),rep(chr+offset+chr.width,2)))
	polygon(x=c(scale.chr,rev(scale.chr)),y=c(rep(chr-offset,2),rep(chr-offset+chr.width,2)))
 })
 
 }
 
 ##next two functions are modified from the plotrix package in R
## http://cran.r-project.org/web/packages/plotrix/index.html
 
 drawSectorAnnulus<-function (angle1, angle2, radius1, radius2, col, angleinc = 0.03,border="lightgrey") 
{
	if(col=="white") return()
	if(abs(angle1-angle2)<.025) border=col
	if(col=="white") border="white"
    if (angle1 > angle2) {
        temp <- angle1
        angle1 <- angle2
        angle2 <- temp
    }
    if (radius1 > radius2) {
        temp <- radius1
        radius1 <- radius2
        radius2 <- temp
    }
    angles <- seq(angle1, angle2, by = abs(angle1-angle2)/30) # angleinc); 
    angles[length(angles)] <- angle2
    xpos <- c(cos(angles) * radius1, cos(rev(angles)) * radius2)
    ypos <- c(sin(angles) * radius1, sin(rev(angles)) * radius2)
    polygon(xpos, ypos, col = col, border = border)
}

my.pie<-function(radial.extents, sector.edges = NULL, sector.colors = NULL,mar = c(2, 2, 3, 2), radial.lim = NULL,clockwise=TRUE,start = 0,add=FALSE,lower.y=-1.5){
    par(mar = mar, pty = "s")
    maxrad <- max(unlist(radial.extents))
 #	recover()
    if (is.null(radial.lim)) 
        radial.lim <- range(radial.extents)
    if (is.null(sector.edges)) {
        if (clockwise) 
            sector.edges <- seq(pi + start, start, length.out = length(radial.extents) + 
                1)
        else sector.edges <- seq(start, 2 * pi + start, length.out = length(radial.extents) + 
            1)
    }
        if(add==FALSE)    plot(0, xlim = c(-maxrad, maxrad), ylim = c(lower.y, 
            maxrad), type = "n", axes = FALSE,xlab="",ylab="")
                nsectors <- length(radial.extents)
if (is.list(radial.extents)) {
        if (is.null(sector.colors)) 
            sector.colors <- rainbow(nsectors)
        for (sector in 1:nsectors) {
            annuli <- radial.extents[[sector]]
            annulus.colors <- sector.colors[[sector]]
            for (annulus in 1:(length(annuli) - 1)) {
                drawSectorAnnulus(sector.edges[[sector]], sector.edges[[sector + 
                  1]], annuli[annulus], annuli[annulus + 1], 
                  annulus.colors[annulus])
            }
        }
    }
}

##overlap<-IBD.overlap.overall(my.blocks.1=family.1.chunks,my.blocks.2=family.2.chunks,meiosis=1)
IBD.overlap.overall<-function(my.blocks.1,my.blocks.2,meiosis,sample.down=FALSE){
	
	num.ancs<-1:2^(meiosis)
	if(sample.down) num.ancs<-sample(num.ancs,20)
	##overlap between 1/2 relatives, e.g. half 1st cousins
	overlap.half<-sapply(num.ancs,function(anc.1){
		sapply(num.ancs,function(anc.2){
			IBD.overlap(relly.pos1=anc.1,relly.pos2=anc.2,my.blocks.1=my.blocks.1,my.blocks.2=my.blocks.2,meiosis=meiosis)
		})
	})

	num.ancs<-1:2^(meiosis)
	mammas<-num.ancs[(num.ancs %% 2)==1]
	pappas<-num.ancs[(num.ancs %% 2)==0]
	par.combos<-1:2^(meiosis-1)
	if(sample.down) par.combos<-sample(par.combos,20)
	##overlap between full relatives, e.g. full cousins
	overlap.full<-sapply(par.combos,function(par.1){
		sapply(par.combos,function(par.2){		
				mat.ibd<-IBD.overlap(relly.pos1=mammas[par.1],relly.pos2=mammas[par.2],my.blocks.1=my.blocks.1,my.blocks.2=my.blocks.2,meiosis=meiosis)
				pat.ibd<-IBD.overlap(relly.pos1=pappas[par.1],relly.pos2=pappas[par.2],my.blocks.1=my.blocks.1,my.blocks.2=my.blocks.2,meiosis=meiosis)
	#			cat(mat.ibd,pat.ibd,"\n")
				return(mat.ibd+pat.ibd)
		})
	})	
	overlap	<-list()
	overlap[["half"]]<-c(overlap.half)
	overlap[["full"]]<-c(overlap.full)
	
	return(overlap)
}


IBD.overlap<-function(relly.pos1,relly.pos2,my.blocks.1,my.blocks.2,meiosis){
	num.ancs<-1:2^(meiosis+1)
	maternal.chr1<-num.ancs[(num.ancs %% 2)==0][relly.pos1]
	paternal.chr1<-num.ancs[(num.ancs %% 2)==1][relly.pos1]

	maternal.chr2<-num.ancs[(num.ancs %% 2)==0][relly.pos2]
	paternal.chr2<-num.ancs[(num.ancs %% 2)==1][relly.pos2]
		
	chr.overlap<-sapply(1:22,function(chr){  
		blocks.1<-my.blocks.1[[meiosis+1]][[maternal.chr1]][[chr]]
		blocks.2<-my.blocks.2[[meiosis+1]][[maternal.chr1]][[chr]]
		num.1<-nrow(blocks.1);num.2<-nrow(blocks.2);
		if(is.null(num.1) | is.null(num.2)){ 
			num.blocks.overlap.1<-0 
			}else{
			heads.tails<-cbind(c(blocks.1$l,blocks.1$r,blocks.2$l,blocks.2$r),c(rep(1,num.1),rep(-1,num.1),rep(1,num.2),rep(-1,num.2)))
			num.blocks.overlap.1<-sum(cumsum(heads.tails[order(heads.tails[,1]),2])==2)
		}
		blocks.1<-my.blocks.1[[meiosis+1]][[paternal.chr2]][[chr]]
		blocks.2<-my.blocks.2[[meiosis+1]][[paternal.chr2]][[chr]]
		num.1<-nrow(blocks.1);num.2<-nrow(blocks.2);
		if(is.null(num.1) | is.null(num.2)){ 
			num.blocks.overlap.2<-0 
			}else{
			heads.tails<-cbind(c(blocks.1$l,blocks.1$r,blocks.2$l,blocks.2$r),c(rep(1,num.1),rep(-1,num.1),rep(1,num.2),rep(-1,num.2)))	
			num.blocks.overlap.2<-sum(cumsum(heads.tails[order(heads.tails[,1]),2])==2)
		}

		return(num.blocks.overlap.1+num.blocks.overlap.2)
	})
#	recover()
	return(sum(chr.overlap))
}

