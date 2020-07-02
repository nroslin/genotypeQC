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

pca<- read.table(input.pca, header=F )
names( pca )[1:5]<-c("FID","IID", "PC1","PC2","PC3" )
pca$ID<- paste( pca$FID,pca$IID, sep=":" )

pca<- merge( ped, pca , by="ID" ) 



# logical indicates which are the samples (ie not 1kg)
samples<- is.na( pca$Continent )



# outliers wrt samples only 



outliers<- 
  apply( pca[ samples ,c( "PC1","PC2","PC3" ) ] , 2, 
      function(v){ 
        bp<- boxplot(v, plot=F, range=3 )
        is.element( v, bp$out )
      }
  )



outliers<- apply( outliers, 1, any ) 


write.table( pca[ samples, c("FID","IID")][outliers,] , paste( output.prefix, ".outliers.txt" , sep="" ), quote=F, col=T, row=F )


col<- as.integer( pca$Continent ) +1 
# black
col[samples]<- 1

pch<- rep( 21, nrow(pca) )
# samples are bullets
pch[ is.na( pca$Continent ) ]<-19 
# outliers are circles with X in them  
pch[samples][ outliers ]<- 13

pdf( paste( output.prefix, ".pcaplot.pdf",sep=""),  width=12, height=18)
#par( mfrow=c(3,1) )
layout( matrix( 1:6, nrow=3, byr=F ) ) 

plot( pca[,c("PC1","PC2") ] , col=col, xlab="PC1", ylab="PC2", pch=pch  ) 
points(  pca[samples ,c("PC1","PC2") ] , col=col[samples], pch=pch[samples]   ) 
legend( x="bottomleft", legend=c("samples", lev<-levels(pca$Continent) ), pch=19, col=1:(1+length( lev )) )
plot( pca[,c("PC1","PC3")] , col=col , xlab="PC1", ylab="PC3", pch=pch  ) 
points(  pca[samples ,c("PC1","PC3")] , col=col[samples], pch=pch[samples]  ) 
plot( pca[,c("PC2","PC3")] , col=col , xlab="PC2", ylab="PC3", pch=pch ) 
points(  pca[samples ,c("PC2","PC3") ] , col=col[samples], pch=pch[samples]  ) 

plot( pca[,c("PC1","PC2") ] , col=col, xlab="PC1", ylab="PC2", pch=pch  ) 
points(  pca[!samples ,c("PC1","PC2") ] , col=col[!samples], pch=pch[!samples]   ) 
legend( x="bottomleft", legend=c("samples", lev<-levels(pca$Continent) ), pch=19, col=1:(1+length( lev )) )
plot( pca[,c("PC1","PC3")] , col=col , xlab="PC1", ylab="PC3", pch=pch  ) 
points(  pca[!samples ,c("PC1","PC3")] , col=col[!samples], pch=pch[!samples]  ) 
plot( pca[,c("PC2","PC3")] , col=col , xlab="PC2", ylab="PC3",pch=pch ) 
points(  pca[!samples ,c("PC2","PC3") ] , col=col[!samples], pch=pch[!samples]  ) 


dev.off()



###############
# distance between samples and 1kg samples. 

sel<-  is.na( pca$Continent )
di<- as.matrix( dist( pca[,-(1:4)] ) )
di<- di[sel, !sel]

# take closest 1kg sample
anc <- pca$Continent[ !sel ]
closest<- apply( di, 1, which.min ) 
sample.anc<- anc[ closest ]

write.table( data.frame( FID=pca$FID[sel], IID=pca$IID[sel], closestAncestry= sample.anc ), 
  paste0(output.prefix,".closestAncestry.txt" ), col=T, row=F, quote=F)












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