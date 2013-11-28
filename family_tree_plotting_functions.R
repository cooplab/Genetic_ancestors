plot.blocks<-function(this.anc,meiosis,chr=1,which.side){	
	blocks<-family.chunks[[meiosis]][[this.anc]][[chr]];
	if(length(blocks)==0){ return()}

	if(which.side=="left") anc.pos<- -1*meiosis-0.5
	if(which.side=="right")  anc.pos<- meiosis + 0.5 
blocks<-(blocks/chr.length)*(-5)-1
#	segments(anc.pos,blocks$l,anc.pos,blocks$r,col=anc.cols[this.anc],lwd=3,lend=2)
	apply(blocks,1,function(x){x= polygon(c(rep(anc.pos,2),rep(anc.pos+.4,2)) , y=c(x,rev(x)),col=anc.cols[this.anc],border="lightgrey")})
}


##draw a chr.

chr.chunks<-function(my.blocks,my.col){
	sapply(1:22,function(chr){  
	   	blocks<-my.blocks[[1]][[chr]]/chr.lengths[1] 
		if(!is.null(nrow(blocks))) apply(blocks,1,function(x){x= polygon(x=c(x,rev(x)),c(rep(chr+offset,2),rep(chr+offset+chr.width,2)),col=my.col)})
		 blocks<-my.blocks[[2]][[chr]]/chr.lengths[1]
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