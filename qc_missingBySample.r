# Created 31 May 2023.
# Last modified:  31 May 2023

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
dev.off()

sessionInfo()
