# Created 20 January 2022.
# Last modified:  21 Jan 2022

# Make some plots and do some sex inference.

infile<-commandArgs(trailingOnly=T)
outfile<-paste0(infile, "_inferredSex.txt")

x<-read.table(paste0(infile, "_sexCheck.txt"), header=T)

#make some plots
pdf(paste0(infile, "_sexCheck.pdf"), height=12, width=18)
par(mfcol=c(2,3))

#plot(x$ChrXhetRate, x$ChrYcr, main="All")
pu<-x[x$PedSex==0,]
pm<-x[x$PedSex==1,]
pf<-x[x$PedSex==2,]

if ( nrow(pu) > 0 ) {
  plot(sort(pu$ChrXhetRate), main="Ped unknown, chrX het", ylim=c(0,0.5))
  plot(sort(pu$ChrYcr), main="Ped unknown, chrY CR", ylim=c(0,1))
} else {
  plot.new()
  plot.new()
}

if ( nrow(pm) > 0 ) {
  plot(sort(pm$ChrXhetRate), main="Ped male, chrX het", ylim=c(0,0.5))
  plot(sort(pm$ChrYcr), main="Ped male, chrY CR", ylim=c(0,1))
} else {
   plot.new()
   plot.new()
}

if ( nrow(pf) > 0 ) {
  plot(sort(pf$ChrXhetRate), main="Ped female, chrX het", ylim=c(0,0.5))
  plot(sort(pf$ChrYcr), main="Ped female, chrY CR", ylim=c(0,1))
} else {
  plot.new()
   plot.new()
}

dev.off()

#sex inference based on X het
ChrXsex<-rep(0, nrow(x))
ChrYsex<-rep(0, nrow(x))
y<-cbind(x, ChrXsex, ChrYsex)
y<-y[order(y$ChrXhetRate),]
d<-diff(y$ChrXhetRate)
w<-which.max(d)  #find the location of the biggest gap

y$ChrXsex[1:w]<-1  #samples with low chrXhet are male
y$ChrXsex[(w+1):nrow(y)]<-2   #samples with high chrXhet are female


#sex inference based on Y cr
y<-y[order(y$ChrYcr),]
d<-diff(y$ChrYcr)
w<-which.max(d)
y$ChrYsex[1:w]<-2   #samples with low chrYcr are female
y$ChrYsex[(w+1):nrow(y)]<-1

#make overall inference
InferredSex<-y$ChrXsex
sel<-y$ChrYsex != y$ChrXsex   #chrX and Y inference conflict
InferredSex[sel]<-0   #set these to missing

sel<-y$PedSex != 0 & y$ChrXsex != y$PedSex   #pedsex differs from chr x sex
InferredSex[sel]<-0   #set these to missing

write.table(cbind(y$FID, y$IID, y$PedSex, InferredSex, y$ChrXsex, y$ChrYsex),
  outfile, row.names=F, col.names=c("FID", "IID", "PedSex",
  "InferredSex", "ChrXsex", "ChrYsex"), quote=F, sep="\t")
