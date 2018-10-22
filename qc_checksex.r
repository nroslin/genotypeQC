FILE<- commandArgs( TRUE )

s<-read.table(FILE, head=T )

s<-s[ order( s$F ),] 
pdf( paste( FILE, ".pdf", sep=""), height=6, width=18 )
par( mfrow=c(1,3) )

if( sum(  s$PEDSEX==0  ) > 0 ){
 plot(  s$F[ s$PEDSEX==0 ]   , col=as.integer(s$STATUS[ s$PEDSEX==0 ] ) , 
        main="unknown sex" , ylim=c(-1,1)  )
} else {
 plot.new()
}
if( sum(  s$PEDSEX==1  ) > 0 ){
 plot(  s$F[ s$PEDSEX==1 ]   , col=as.integer(s$STATUS[ s$PEDSEX==1 ] ) , 
        main="males" , ylim=c(-1,1)  )
} else {
 plot.new()
}

if( sum(  s$PEDSEX==2  ) > 0 ){
 plot(  s$F[ s$PEDSEX==2 ]   , col=as.integer(s$STATUS[ s$PEDSEX==2 ] ) , 
        main="females", ylim=c(-1,1)  )
} else {
 plot.new()
}

