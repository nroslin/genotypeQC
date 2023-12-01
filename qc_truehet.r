ARGS<- commandArgs( TRUE )

prefix<- ARGS[1]
bprange <- as.integer( ARGS[2] )

s<-read.table( paste0(prefix, "_autoHet.txt"), head=T )

pdf( paste( prefix , "_autoHet.pdf", sep=""), height=6, width=18 )
par( mfrow=c(1,3) )

s$AutoHet<-s$Nhet / s$Ngtypes
bp<- boxplot( s$AutoHet, range=bprange , plot=F  )
lo<-  bp$stats[1,1]
hi<-  bp$stats[5,1] 
hist(s$AutoHet, main="", xlab="Autosomal heterozygosity", col="grey")
plot(sort(s$AutoHet), ylab="Autosomal heterozygosity")
abline( h= bp$stats[c(1,5),1], lty=3 )
boxplot( s$AutoHet, range=bprange )


sel<- s$AutoHet > hi | s$AutoHet < lo 
if( sum( sel ) > 0 ){
 df<- cbind( s[sel, 2:3 ], SOURCE="HET" )
 write.table( df, paste0(prefix, ".remove.txt"), quote=F, col=F, row=F, append=T )
}

mtext(paste(date()), outer=T, side=3, adj=1, line=-1.2, cex=0.7)
mtext(paste0(system("pwd", intern=T), "/", prefix, "_autoHet.pdf"), outer=T,
  side=1, adj=0, line=-1, cex=0.7, font=3)

dev.off()


sessionInfo()








#############################################################################
#   Copyright 2019 Mathieu Lemire
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#############################################################################
