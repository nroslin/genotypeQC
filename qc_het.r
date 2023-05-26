ARGS<- commandArgs( TRUE )

input<- ARGS[1]
output <-ARGS[2]
bprange <- as.integer( ARGS[3] )

s<-read.table( input , head=T )

pdf( paste( input , ".pdf", sep=""), height=6, width=18 )
par( mfrow=c(1,3) )

bp<- boxplot( s$F, range=bprange , plot=F  )
lo<-  bp$stats[1,1]
hi<-  bp$stats[5,1] 
plot(  sort( s$F ) , main="het/F statistic" )
abline( h= bp$stats[c(1,5),1], lty=3 )
boxplot( s$F, range=bprange )

hist(s$F, col="grey", main="", xlab="het/F-statistic")

sel<- s$F > hi | s$F < lo 
if( sum( sel ) > 0 ){
 df<- cbind( s[sel, 1:2 ], SOURCE="HET" )
 write.table( df, output, quote=F, col=F, row=F, append=T )
}

dev.off()











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
