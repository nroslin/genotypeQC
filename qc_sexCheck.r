# Created 20 January 2022.
# Last modified:  25 May 2023

# Make some plots and do some sex inference.
# Estimate how many of each sex chromosome a person carries (0, 1, 2, 
# undetermined), then use this to estimate sex (male, female, unclear).  The
# ones that are unclear can be examined for possible anomalies (eg. XXY).
# Use fixed threshold to establish counts.

infile<-commandArgs(trailingOnly=T)
outfile<-paste0(infile, "_inferredSex.txt")


x<-read.table(paste0(infile, "_sexCheck.txt"), header=T, as.is=T)

#make some plots
pdf(paste0(infile, "_sexCheck.pdf"), height=5, width=15)
par(mfrow=c(1,3))
#chr X het vs. chr Y call rate
plot(x$ChrXhetRate, x$ChrYcr, main=paste(infile, ": All samples"), pch=x$PedSex+1, 
  col=x$PedSex+1, xlab="ChrX heterozygosity", ylab="ChrY call rate")
segments(0.10, 0.4, 0.3, 0.4, lty=2)
segments(0.10, 0.4, 0.10, -1, lty=2)
segments(0.3, 0.4, 0.3, -1, lty=2)
segments(-1, 0.9, 0.05, 0.9, lty=2)
segments(0.05, 0.9, 0.05, 1.1, lty=2)

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

#chrY
ny<-rep(NA, nrow(x))
sel<-x$ChrYcr>0.9   #one Y chr if chrY(cr)>0.9
ny[sel]<-1
sel<-x$ChrYcr<0.4   #no Y chr if chrY(cr)<0.4  (chip-dependent)
ny[sel]<-0
#rest are NA


### inferred sex ###
inf.sex<-rep(NA, nrow(x))
sel<-nx==1 & ny==1   #1 X, 1 Y (male)
inf.sex[sel]<-1
sel<-nx==2 & ny==0   #2 X, 0 Y (female)
inf.sex[sel]<-2


#write it out, without doing anything about conflict with input sex/gender
write.table(cbind(x$FID, x$IID, x$PedSex, inf.sex, nx, ny), outfile, row.names=F,
  col.names=c("FID", "IID", "PedSex", "InferredSex", "NchrX", "NchrY"),
  quote=F, sep="\t")

#par(mfcol=c(2,3))

#pu<-x[x$PedSex==0,]
#pm<-x[x$PedSex==1,]
#pf<-x[x$PedSex==2,]

#if ( nrow(pu) > 0 ) {
#  plot(sort(pu$ChrXhetRate), main="Ped unknown, chrX het", ylim=c(0,0.5))
#  plot(sort(pu$ChrYcr), main="Ped unknown, chrY CR", ylim=c(0,1))
#} else {
#  plot.new()
#  plot.new()
#}

#if ( nrow(pm) > 0 ) {
#  plot(sort(pm$ChrXhetRate), main="Ped male, chrX het", ylim=c(0,0.5))
#  plot(sort(pm$ChrYcr), main="Ped male, chrY CR", ylim=c(0,1))
#} else {
#   plot.new()
#   plot.new()
#}

#if ( nrow(pf) > 0 ) {
#  plot(sort(pf$ChrXhetRate), main="Ped female, chrX het", ylim=c(0,0.5))
#  plot(sort(pf$ChrYcr), main="Ped female, chrY CR", ylim=c(0,1))
#} else {
#  plot.new()
#   plot.new()
#}

#dev.off()

#sex inference based on X het
#ChrXsex<-rep(0, nrow(x))
#ChrYsex<-rep(0, nrow(x))
#y<-cbind(x, ChrXsex, ChrYsex)
#y<-y[order(y$ChrXhetRate),]
#d<-diff(y$ChrXhetRate)
#w<-which.max(d)  #find the location of the biggest gap

#y$ChrXsex[1:w]<-1  #samples with low chrXhet are male
#y$ChrXsex[(w+1):nrow(y)]<-2   #samples with high chrXhet are female


#sex inference based on Y cr
#y<-y[order(y$ChrYcr),]
#d<-diff(y$ChrYcr)
#w<-which.max(d)
#y$ChrYsex[1:w]<-2   #samples with low chrYcr are female
#y$ChrYsex[(w+1):nrow(y)]<-1

#make overall inference
#InferredSex<-y$ChrXsex
#sel<-y$ChrYsex != y$ChrXsex   #chrX and Y inference conflict
#InferredSex[sel]<-0   #set these to missing

#sel<-y$PedSex != 0 & y$ChrXsex != y$PedSex   #pedsex differs from chr x sex
#InferredSex[sel]<-0   #set these to missing

#write.table(cbind(y$FID, y$IID, y$PedSex, InferredSex, y$ChrXsex, y$ChrYsex),
#  outfile, row.names=F, col.names=c("FID", "IID", "PedSex",
#  "InferredSex", "ChrXsex", "ChrYsex"), quote=F, sep="\t")
