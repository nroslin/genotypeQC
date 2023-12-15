# Created 31 May 2023.
# Last modified:  13 Dec 2023

# Draw qq plots for batch effects analysis.
# Input file is a .model file from PLINK with standard headers

ARGS<-commandArgs(T)

infile<-ARGS[1]
label<-ARGS[2]   #batch name
infile
label

source("/hpf/projects/arnold/users/nroslin/Scripts/R/qqman.r")
options(bitmapType="cairo")
png(paste0("qq_", label, ".png"), width=800, height=800)
#par(mfrow=c(1,2))

x<-read.table(infile, header=T)
y<-x[x$TEST=="ADD",]  #just additive genetic results
#manhattan(x)
qq(y$P, cex.axis=1.5, pt.cex=0.7, cex.lab=1.5, cex.main=1.5)
mtext(paste(label), adj=0, font=2, side=3)
dev.off()

