source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

D1=read.table('MGH36_mat.txt',header=T,row.names=1,sep='\t')
D2=read.table('MGH53_mat.txt',header=T,row.names=1,sep='\t')
D3=read.table('MGH54_mat.txt',header=T,row.names=1,sep='\t')
D4=read.table('MGH60_mat.txt',header=T,row.names=1,sep='\t')
D5=read.table('MGH93_mat.txt',header=T,row.names=1,sep='\t')
D6=read.table('MGH97_mat.txt',header=T,row.names=1,sep='\t')

saveRDS(D1,file='MGH36.RDS')
saveRDS(D2,file='MGH53.RDS')
saveRDS(D3,file='MGH54.RDS')
saveRDS(D4,file='MGH60.RDS')
saveRDS(D5,file='MGH93.RDS')
saveRDS(D6,file='MGH97.RDS')



BATCH=c(rep('D1',ncol(D1)),
      rep('D2',ncol(D2)),
      rep('D3',ncol(D3)),
      rep('D4',ncol(D4)),
      rep('D5',ncol(D5)),
      rep('D6',ncol(D6)) )


D12=.simple_combine(D1,D2)$combine
D34=.simple_combine(D3,D4)$combine
D56=.simple_combine(D5,D6)$combine
D1234=.simple_combine(D12,D34)$combine
D123456=.simple_combine(D1234,D56)$combine

DATA=D123456

dim(DATA) #23686  4347
rm(D1)
rm(D2)
rm(D3)
rm(D4)
rm(D5)
rm(D6)
rm(D12)
rm(D34)
rm(D56)
rm(D1234)
rm(D123456)



mybeer=MBEER(DATA, BATCH, CNUM=10, PCNUM=50,CPU=4)

pbmc=mybeer$seurat

ALLPC=c(1:20)
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = ALLPC, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)


PCUSE <- which(mybeer$cor > min(0.7, median(mybeer$cor))  & mybeer$fdr<0.05 & c(1:50)<=20  )
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims.use = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)


