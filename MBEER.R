source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')





BEER <- function(DATA, TAG, CNUM=10, PCNUM=50, VPCOR=0, CPU=4, print_step=10){
  
  
    
    RESULT=list()
    library(Seurat)
    source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BEER start!')
    print(Sys.time())
    D1=D1
    D2=D2
    CNUM=CNUM
    PCNUM=PCNUM
    print_step=print_step
    
    print('############################################################################')
    print('MainStep1.Combine Data...')
    print('############################################################################')
    EXP=.simple_combine(D1,D2)$combine
    pbmc=CreateSeuratObject(raw.data = EXP, min.cells = 0, min.genes = 0, project = "ALL") 
    
    print('############################################################################')
    print('MainStep2.Preprocess Data...')
    print('############################################################################')
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc <- FindVariableGenes(object = pbmc, do.plot=F,mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    #length(x = pbmc@var.genes)
    pbmc <- ScaleData(object = pbmc, genes.use=pbmc@var.genes, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
    pbmc <- RunPCA(object = pbmc, pcs.compute=PCNUM,pc.genes = pbmc@var.genes, do.print =F)
    
    print('############################################################################')
    print('MainStep3.Convert to one-dimension...')
    print('############################################################################')
    D1X=.data2one(D1, pbmc@var.genes, CPU, PCNUM)
    D2X=.data2one(D2, pbmc@var.genes, CPU, PCNUM)
    G1=.getGroup(D1X,'D1',CNUM)
    G2=.getGroup(D2X,'D2',CNUM)
    GROUP=c(G1,G2)
    BATCH=c(rep('D1',ncol(D1)),rep('D2',ncol(D2)))
    pbmc@meta.data$group=GROUP
    pbmc@meta.data$batch=BATCH
    
    
    print('############################################################################')
    print('MainStep4.Get Valid Pairs...')
    print('############################################################################')
    VP_OUT=.getValidpair(D1, G1, D2, G2, CPU, method='kendall', print_step)
    #VP_OUT=.getValidpair(D1, G1, D2, G2, 4, 'kendall', 10)
    VP=VP_OUT$vp
    ##########################
    NROW_VP=nrow(VP)
    print('n(Validpair):')
    print(NROW_VP)
    #if(NROW_VP<=1 | is.null(NROW_VP) ){print('Please set a smaller CNUM !!!')}
    if(NROW_VP<=1 | is.null(NROW_VP) ){return(message("Please set a smaller CNUM !!!"))}
    ##########################
    VP=VP[which(VP_OUT$cor>=VPCOR),]
    MAP=rep('NA',length(GROUP))
    MAP[which(GROUP %in% VP[,1])]='D1'
    MAP[which(GROUP %in% VP[,2])]='D2'
    pbmc@meta.data$map=MAP
    
    print('############################################################################')
    print('MainStep5.Detect subspaces with batch effect...')
    print('############################################################################')
    DR=pbmc@dr$pca@cell.embeddings 
    B1index=which(BATCH=='D1')
    B2index=which(BATCH=='D2')
    OUT=.evaluateBatcheffect(DR, B1index, B2index, GROUP, VP)
    
    ########################## 
    RESULT$seurat=pbmc
    RESULT$vp=VP
    RESULT$vpcor=VP_OUT$cor
    RESULT$d1x=D1X
    RESULT$d2x=D2X
    RESULT$g1=G1
    RESULT$g2=G2
    RESULT$cor=OUT$cor
    RESULT$pv=OUT$pv
    RESULT$fdr=OUT$fdr
    
    #RESULT$pcuse=PCUSE
    print('############################################################################')
    print('BEER cheers !!! All main steps finished.')
    print('############################################################################')
    print(Sys.time())
    return(RESULT)
    }
