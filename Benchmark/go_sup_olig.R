MGH36=read.table('MGH36.txt',header=T,row.names=1,sep='\t')
MGH53=read.table('MGH53.txt',header=T,row.names=1,sep='\t')
MGH54=read.table('MGH54.txt',header=T,row.names=1,sep='\t')
MGH60=read.table('MGH60.txt',header=T,row.names=1,sep='\t')
MGH93=read.table('MGH93.txt',header=T,row.names=1,sep='\t')
MGH97=read.table('MGH97.txt',header=T,row.names=1,sep='\t')


library(Seurat)
source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')


D1 <- MGH36
D2 <- MGH53
mybeer <- BEER(D1, D2, CNUM=10, PCNUM=50, CPU=2)










