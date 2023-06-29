# Created 20 January 2022.
# Last modified:  29 Jun 2023

# Make some plots and do some sex inference.
# Estimate how many of each sex chromosome a person carries (0, 1, 2, 
# undetermined), then use this to estimate sex (male, female, unclear).  The
# ones that are unclear can be examined for possible anomalies (eg. XXY).
# Use fixed threshold to establish counts.

infile<-commandArgs(trailingOnly=T)
outfile<-paste0(infile, "_inferredSex.txt")


x<-read.table(paste0(infile, "_sexCheck.txt"), header=T, as.is=T)




### estimate counts ###
#chrX
nx<-rep(NA, nrow(x))
sel<-x$ChrXhetRate>0.1 & x$ChrXhetRate<0.3  #two X chr if 0.1<chrX(het)<0.3
nx[sel]<-2
sel<-x$ChrXhetRate<0.05   #one X chr if chrX(het)<0.05
nx[sel]<-1
sel<-x$ChrXhetRate>0.3   #three X chr if chrX(het)>0.3  ###adjust threshold?
nx[sel]<-3
#rest are NA
x$NchrX<-nx

#chrY
ny<-rep(NA, nrow(x))
sel<-x$ChrYcr>0.75   #one Y chr if chrY(cr)>0.75
ny[sel]<-1
sel<-x$ChrYcr<0.4   #no Y chr if chrY(cr)<0.4  (chip-dependent)
ny[sel]<-0
#rest are NA
x$NchrY<-ny


### inferred sex ###
inf.sex<-rep(NA, nrow(x))
sel<-nx==1 & ny==1   #1 X, 1 Y (male)
inf.sex[sel]<-1
sel<-nx==2 & ny==0   #2 X, 0 Y (female)
inf.sex[sel]<-2
x$InferredSex<-inf.sex



#write it out, without doing anything about conflict with input sex/gender
write.table(x[,c("FID", "IID", "PedSex", "InferredSex", "NchrX", "NchrY")], 
  outfile, row.names=F, quote=F, sep="\t")


#make some plots
if ( sum(!is.na(x$ChrXcr)) > 0 ) {  #chrX call rate has been calculated;
					#means extra  plots will be drawn
  pdf(paste0(infile, "_sexCheck.pdf"), height=5, width=20)
  par(mfrow=c(1,4))
} else {
  pdf(paste0(infile, "_sexCheck.pdf"), height=5, width=10)
  par(mfrow=c(1,2))
}


#chr X het vs. chr Y call rate
plot(x$ChrXhetRate, x$ChrYcr, main=paste0(infile, ": ", nrow(x), " samples"), 
  pch=x$PedSex+1, 
  col=x$PedSex+1, xlab="ChrX heterozygosity", ylab="ChrY call rate")
segments(0.10, 0.4, 0.3, 0.4, lty=2)
segments(0.10, 0.4, 0.10, -1, lty=2)
segments(0.3, 0.4, 0.3, -1, lty=2)
segments(-1, 0.75, 0.05, 0.9, lty=2)
segments(0.05, 0.75, 0.05, 1.1, lty=2)

#repeat for problematic samples and ones where pedsex is unknown (coded as 0)
problem<-x$PedSex==0 | x$PedSex != x$InferredSex | is.na(x$InferredSex)
plot(x[problem, "ChrXhetRate"], x[problem,"ChrYcr"], main=paste0(infile, ": ",
  sum(problem), " atypical samples"), pch=x[problem,"PedSex"]+1, 
  col=x[problem,"PedSex"]+1, xlab="ChrX heterozygosity", ylab="ChrY call rate")
segments(0.10, 0.4, 0.3, 0.4, lty=2)
segments(0.10, 0.4, 0.10, -1, lty=2)
segments(0.3, 0.4, 0.3, -1, lty=2)
segments(-1, 0.75, 0.05, 0.9, lty=2)
segments(0.05, 0.75, 0.05, 1.1, lty=2)

### two optional plots
#chrX het vs. chr X call rate
#chrX call rate may not have been calculated, so check if all missing first
if ( sum(!is.na(x$ChrXcr)) > 0 ) {
plot(x$ChrXhetRate, x$ChrXcr, main=paste(infile, ": All samples"), pch=x$PedSex+1, 
  col=x$PedSex+1, xlab="ChrX heterozygosity", ylab="ChrX call rate")

#chrX call rate vs. chrY call rate
plot(x$ChrXcr, x$ChrYcr, main=paste(infile, ": All samples"), pch=x$PedSex+1, 
  col=x$PedSex+1, xlab="ChrX call rate", ylab="ChrY call rate")
}

dev.off()

