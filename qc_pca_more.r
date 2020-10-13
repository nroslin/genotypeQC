
Npc<- 10 


ARGS<- commandArgs( TRUE )

# input.prefix is a raw file produced with --recodeA option of plink 
input.pca <- ARGS[1]
input.fam<- ARGS[2] 
output.prefix  <-ARGS[3]

kgsamples <- read.table("/hpf/projects/arnold/references/www.tcag.ca/documents/tools/sampleTable.txt", header=T )
pop<- kgsamples[, c("ID","Continent")]
names(pop)[1]<-"IID"

gen <-read.table( input.fam , head=F )
ped<- data.frame( gen[,1:2], o=1:nrow( gen ) )
names( ped )[1:2]<-c("FID","IID" )
ped<- merge( ped, pop, all.x=T , by="IID"  ) 
ped<- ped[ order(ped$o), -3 ] 
ped$ID<- paste( ped$FID,ped$IID, sep=":" )
ped<- ped[ , -(1:2)]

input.val<- sub( "eigenvec","eigenval", input.pca )
val<- read.table(input.val )

pca<- read.table(input.pca, header=F )
names( pca )[1:2]<-c("FID","IID" )
names( pca )[ 3:ncol(pca) ]<- paste0( "PC", 3:ncol(pca) -2 ) 

if( Npc > ncol(pca) -2 ){ Npc= ncol(pca) -2 }

pca$ID<- paste( pca$FID,pca$IID, sep=":" )



pca<- merge( ped, pca , by="ID" ) 




# logical indicates which are the samples (ie not 1kg)
samples<- is.na( pca$Continent )



# outliers wrt samples only 

col<- as.integer( pca$Continent ) +1 
# black
col[samples]<- 1

pch<- rep( 21, nrow(pca) )
# samples are bullets
pch[ is.na( pca$Continent ) ]<-19 

cex<- rep( 1,  nrow(pca) )
cex[ is.na( pca$Continent ) ]<- .2 


pcnames<- paste0( "PC", 1:Npc   )

for( i in 1:( length( pcnames )-1 ) ){
 for( j in (i+1):length( pcnames ) ){

cat(i," ",j,"\n") 
pdf( paste( output.prefix,".",i,".",j, ".pcaplot.pdf",sep=""),  width=12, height=6)

layout( matrix( 1:2, nrow=1, byr=F ) ) 
plot( pca[, pcnames[ c( i,j ) ] ] , col=col, xlab=pcnames[i], ylab=pcnames[j] , pch=pch , cex=cex  ) 
points(  pca[samples , pcnames[ c( i,j ) ] ]  , col=col[samples], pch=pch[samples] , cex=cex[samples]  ) 
legend( x="bottomleft", legend=c("samples", lev<-levels(pca$Continent) ), pch=19, col=1:(1+length( lev )) )
plot( pca[, pcnames[ c( i,j ) ] ] , col=col, xlab=pcnames[i], ylab=pcnames[j] , pch=pch , cex=cex  )
points(  pca[!samples , pcnames[ c( i,j ) ] ]  , col=col[!samples], pch=pch[! samples] , cex=cex[! samples]  )
legend( x="topleft", legend=c("samples", lev<-levels(pca$Continent) ), pch=19, col=1:(1+length( lev )) )

dev.off()

}
} 















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