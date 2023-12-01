# Created 31 May 2023.
# Last modified:  01 Dec 2023

# Plot call rate per sample

ARGS<-commandArgs(T)

prefix<-ARGS[1]  #input file prefix
thresh<-ARGS[2]   #missingness threshold

thresh<-as.numeric(thresh)

x<-read.table(paste0(prefix, "_missing_step2.imiss"), header=T, as.is=T)
cr<-1-x$F_MISS   #call rate
crrange<-c(min(c(cr, 1-thresh)), 1)

pdf(paste0(prefix, "_sampleCallRate.pdf"), height=6, width=6)
  plot((1:nrow(x)-1)/nrow(x), sort(cr), ylab="Call Rate", ylim=crrange,
	xlab="Quantile", main=paste0(prefix, ": ", nrow(x), " samples"))
  grid()
  abline(h=1-thresh, col="blue", lty=2)
mtext(paste(date()), outer=T, side=3, adj=1, line=-1.2, cex=0.7)
mtext(paste0(system("pwd", intern=T), "/", prefix, "_sampleCallRate.pdf"), 
  outer=T, side=1, adj=0, line=-1, cex=0.7, font=3)
dev.off()

sessionInfo()
