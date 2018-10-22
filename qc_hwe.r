
ARGS<- commandArgs( TRUE )

input<- ARGS[1]
output <-ARGS[2]

a<-read.table(input , head=T )

a<- a[ a$TEST=="ALL",]

fdr<- p.adjust( a$P, method="fdr" )

sel<- fdr < 0.01 
label<- paste( "HWE:FDR-LOG10=",-floor( log10( fdr[sel] ) +.5 ),sep="" )
write.table(data.frame( SNP= a$SNP[sel], SOURCE=label), output,
  quote=F, col=F, row=F, append=T )


 




