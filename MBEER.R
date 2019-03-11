source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')
source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')

D1=read.table('MGH36_mat.txt',header=T,row.names=1,sep='\t')
D2=read.table('MGH53_mat.txt',header=T,row.names=1,sep='\t')
D3=read.table('MGH54_mat.txt',header=T,row.names=1,sep='\t')
D4=read.table('MGH60_mat.txt',header=T,row.names=1,sep='\t')
D5=read.table('MGH93_mat.txt',header=T,row.names=1,sep='\t')
D6=read.table('MGH97_mat.txt',header=T,row.names=1,sep='\t')

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







MBEER <- function(DATA, BATCH, CNUM=10, PCNUM=50, CPU=4, print_step=10){
  
    RESULT=list()
    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BEER start!')
    print(Sys.time())
    DATA=DATA
    BATCH=BATCH
    CNUM=CNUM
    PCNUM=PCNUM
    print_step=print_step
    
    print('############################################################################')
    print('MainStep1.Preprocess...')
    print('############################################################################')
    TABLE=table(BATCH)  
    MAXBATCH=rownames(TABLE)[which(TABLE==max(TABLE))[1]]
    ############################################################################
    ############################################################################
    EXP = DATA
    pbmc=CreateSeuratObject(raw.data = EXP, min.cells = 0, min.genes = 0, project = "ALL") 
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    #length(x = pbmc@var.genes)
    pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
    pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)
    pbmc@meta.data$batch=BATCH
      
    PAIR=c()
    for(one in rownames(TABLE)){
        if(one != MAXBATCH){
            PAIR=cbind(PAIR, c(MAXBATCH, one))
            }
        }
    PAIR=t(PAIR)
      
    print('############################################################################')
    print('MainStep2. Analyze Each Pair of Batches...')
    print('############################################################################')
    print('Total Number of Batch Pairs:')
    print(nrow(PAIR))
    COR=c()
    PV=c()
    FDR=c()
      
      
    MAX_D1=EXP[,which(BATCH == MAXBATCH)]
    MAX_D1X=.data2one(MAX_D1, pbmc@var.genes, CPU, PCNUM)  
    MAX_G1=.getGroup(MAX_D1X,'D1',CNUM)
    DR=pbmc@dr$pca@cell.embeddings 
      
    i=1
    while(i<=nrow(PAIR)){
          
    this_pair=PAIR[i,]
    print(this_pair)
    this_D2=EXP[,which(BATCH == this_pair[2])]
    this_D2X=.data2one(this_D2, pbmc@var.genes, CPU, PCNUM)         
    this_G2=.getGroup(this_D2X,'D2',CNUM)
    this_GROUP=c(MAX_G1, this_G2)
    this_BATCH=c(rep('D1',ncol(MAX_D1)),rep('D2',ncol(this_D2))) 
    this_VP_OUT=.getValidpair(MAX_D1, MAX_G1, this_D2, this_G2, CPU, method='kendall', print_step)  
    #this_VP_OUT=.getValidpair(MAX_D1, MAX_G1, this_D2, this_G2, CPU, method='kendall', 10)
    this_VP=this_VP_OUT$vp
    this_NROW_VP=nrow(this_VP)
    print('n(Validpair):')
    print(this_NROW_VP)
    if(this_NROW_VP<=1 | is.null(this_NROW_VP) ){return(message("Please set a smaller CNUM !!!"))}     
    ##########################
    this_DR=DR[c(which(BATCH %in% this_pair[1]), which(BATCH %in% this_pair[2])),]
    this_B1index=which(this_BATCH=='D1')
    this_B2index=which(this_BATCH=='D2')
    this_OUT=.evaluateBatcheffect(this_DR, this_B1index, this_B2index, this_GROUP, this_VP)      
    #################
    COR=cbind(COR,this_OUT$cor)
    PV=cbind(PV,this_OUT$pv)
    FDR=cbind(FDR, this_OUT$fdr)
    print('Solved Number:')
    print(i)
          
    i=i+1}
    
    print('############################################################################')
    print('MainStep3. Output')
    print('############################################################################')
    
    ########################## 
    RESULT$seurat=pbmc
    RESULT$COR=COR
    RESULT$PV=PV
    RESULT$FDR=OUT$FDR
    RESULT$cor=apply(COR, 1, mean)
    RESULT$pv=apply(PV, 1, mean)
    RESULT$fdr=apply(FDR, 1, mean)
    #RESULT$pcuse=PCUSE
    print('############################################################################')
    print('BEER cheers !!! All main steps finished.')
    print('############################################################################')
    print(Sys.time())
    return(RESULT)
    }


OUT=MBEER(DATA, BATCH, CNUM=10, PCNUM=50,CPU=4)
