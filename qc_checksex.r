FILE<- commandArgs( TRUE )

out.inferred <- sub("_check-sex.sexcheck",".inferredsex.txt", FILE) 

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


# creating a file of infered sex 
d<- diff( s$F )
w<- which.max(d)
sel<- s$SNPSEX[1:w] == 0
s$SNPSEX[1:w][sel]<- 2 


sel<- s$SNPSEX[(w+1):nrow(s)] == 0
s$SNPSEX[(w+1):nrow(s)][sel]<- 1 

write.table(s[,c("FID","IID","PEDSEX","SNPSEX")], out.inferred ,col=T,row=F,quote=F )













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