
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