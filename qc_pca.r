ARGS<- commandArgs( TRUE )

# input is a raw file produced with --recodeA option of plink 
input<- ARGS[1]
output.prefix  <-ARGS[2]

kgsamples <- read.table("/hpf/projects/arnold/references/www.tcag.ca/documents/tools/sampleTable.txt", header=T )
pop<- kgsamples[, c("ID","Continent")]
names(pop)[1]<-"IID"

gen <-read.table( input , head=T )
ped<- data.frame( gen[,1:6], o=1:nrow( gen ) )
ped<- merge( ped, pop, all.x=T  ) 
ped<- ped[ order(ped$o), -7 ] 

gen<- gen[,-(1:6)] 

impute<-function(v){
 sel<- is.na(v)
 v[sel]<- mean( v[!sel] )
 return( v )
}

gen.imputed<- apply( gen, 2, impute ) 

system.time( pr<- prcomp( t( gen.imputed ), scale.=F, center=T  , retx=F ) )
pr.ped<- list( pr, ped )
save( pr.ped , file=paste( output.prefix,".prcomp.rda", sep="") ) 


# logical indicates which are the samples (ie not 1kg)
samples<- is.na( ped$Continent )

# outliers wrt samples only 

outliers<- 
  apply( pr$rotation[ samples , 1:3] , 2, 
      function(v){ 
        bp<- boxplot(v, plot=F, range=3 )
        is.element( v, bp$out )
      }
  )

outliers<- apply( outliers, 1, any ) 

write.table( ped[ samples, c("FID","IID")][outliers,] , paste( output.prefix, ".outliers.txt" , sep="" ), quote=F, col=T, row=F )


col<- as.integer( ped$Continent ) +1 
# black
col[samples]<- 1

pch<- rep( 21, nrow(ped) )
# samples are bullets
pch[ is.na( ped$Continent ) ]<-19 
# outliers are circles with X in them  
pch[samples][ outliers ]<- 13

pdf( paste( output.prefix, ".pcaplot.pdf",sep=""),  width=6, height=18)
par( mfrow=c(3,1) )
plot( pr$rotation[,1:2] , col=col, xlab="PC1", ylab="PC2", pch=pch  ) 
points(  pr$rotation[samples ,1:2] , col=col, pch=pch[samples]   ) 
legend( x="bottomleft", legend=c("samples", lev<-levels(ped$Continent) ), pch=19, col=1:(1+length( lev )) )
plot( pr$rotation[,c(1,3)] , col=col , xlab="PC1", ylab="PC3", pch=pch  ) 
points(  pr$rotation[samples ,c(1,3)] , col=col, pch=pch[samples]  ) 
plot( pr$rotation[,2:3] , col=col , xlab="PC2", ylab="PC3") 
points(  pr$rotation[samples ,2:3] , col=col, pch=pch[samples]  ) 
dev.off()




