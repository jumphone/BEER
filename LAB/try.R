library(Seurat)
source('BEER.R')


D1=readRDS('MGH36.RDS')
D2=readRDS('MGH53.RDS')
D3=readRDS('MGH54.RDS')
D4=readRDS('MGH60.RDS')
D5=readRDS('MGH93.RDS')
D6=readRDS('MGH97.RDS')

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

rm(D1);rm(D2);rm(D3);rm(D4);rm(D5);rm(D6)


mybeer=ProBEER(DATA,BATCH)

PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL)

pbmc <- mybeer$seurat
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)





pbmc=readRDS('RAW.RDS')
source('BEER_NEW.R')
EXP=as.matrix(pbmc@assays$RNA@counts)
BATCH=pbmc@meta.data$batch


mybeer=ProBEER(EXP,BATCH)





pbmc <- mybeer$seurat
plot(mybeer$cor)
plot(mybeer$lcor)
plot(rank(mybeer$cor),rank(mybeer$lcor),pch=16)




PCUSE=which(rank(mybeer$cor)>=length(mybeer$cor)/2 & rank(mybeer$lcor) >=length(mybeer$cor)/2)
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL)



pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)
DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)
#DimPlot(pbmc, reduction.use='umap', group.by='map', pt.size=0.1)










EXP=as.matrix(pbmc@assays$RNA@counts)
BATCH=pbmc@meta.data$batch
GN=2000
MTTAG="^mt-"
PCNUM=50
SEED=123
CNUM=50
UBATCH=unique(BATCH)
print_step=10
method='kendall'





.evaluateNew <- function(DR, GROUP, VP){
    
    OUT=list() 
    VALID_PAIR=VP
    ALL_COR=c()   
    ALL_PV=c() 
 
    print('Start')
    THIS_DR=1
    while(THIS_DR<=ncol(DR)){
        THIS_PC = DR[,THIS_DR]
        
        lst1_quantile=c()
        lst2_quantile=c()
        i=1
        while(i<=nrow(VALID_PAIR)){
            this_pair=VALID_PAIR[i,]
            this_index1=which(GROUP %in% this_pair[1])
            this_index2=which(GROUP %in% this_pair[2])                   
            lst1_quantile=c(lst1_quantile,quantile(DR[this_index1,THIS_DR]))
            lst2_quantile=c(lst2_quantile,quantile(DR[this_index2,THIS_DR]))
                       
            i=i+1}
        
        this_test=cor.test(lst1_quantile, lst2_quantile, method='kendall')
                
        this_cor=this_test$estimate
        this_pv=this_test$p.value
         
        ALL_COR=c(ALL_COR, this_cor)
        ALL_PV=c(ALL_PV, this_pv) 
        print(THIS_DR)
        
        THIS_DR=THIS_DR+1}
    
    
    OUT$cor=ALL_COR
    OUT$pv=ALL_PV
    OUT$fdr=p.adjust(ALL_PV,method='fdr')
    print('Finished!!!')
    return(OUT)
   }







ProBEER <- function(DATA, BATCH, MAXBATCH="", CNUM=50, PCNUM=50, GN=2000, CPU=4, print_step=10, SEED=123,MTTAG="^MT-", REGBATCH=FALSE){

    set.seed( SEED)
    RESULT=list()
    library(Seurat)
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BEER start!')
    print(Sys.time())
    DATA=DATA
    BATCH=BATCH
    CNUM=CNUM
    PCNUM=PCNUM
    MTTAG=MTTAG
    GN=GN
    print_step=print_step
    
    VARG=c()
    for(this_batch in UBATCH){
        this_pbmc=CreateSeuratObject(counts = EXP[,which(BATCH==this_batch)], min.cells = 0, 
                                 min.features = 0, project = this_batch)
        this_pbmc <- NormalizeData(object = this_pbmc, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
        this_pbmc <- FindVariableFeatures(object = this_pbmc, selection.method = "vst", nfeatures = GN)  
        this_varg=VariableFeatures(object = this_pbmc)
        VARG=c(VARG, this_varg)
        }
    VARG=unique(VARG)



    pbmc=CreateSeuratObject(counts = EXP, min.cells = 0, min.features = 0, project = "ALL") 
    pbmc@meta.data$batch=BATCH
    VariableFeatures(object = pbmc)=VARG


    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = MTTAG)

    CPU=4
    pbmc <- ScaleData(object = pbmc, features = VARG, vars.to.regress = c("nCount_RNA","percent.mt"), num.cores=CPU, do.par=TRUE)
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=1)


    ONE=pbmc@reductions$umap@cell.embeddings[,1]
    GROUP=rep('NA',length(BATCH))
    for(this_batch in UBATCH){
        this_index=which(BATCH==this_batch)
        this_one=ONE[this_index]
        this_group=.getGroup(this_one,this_batch,CNUM)
        GROUP[this_index]=this_group
    }

    pbmc@meta.data$group=GROUP
    VP=c()
    i=1
    while(i<length(UBATCH)){
        j=i+1
        batch1=UBATCH[i]
        batch1_index=which(BATCH==batch1)
        exp1=as.matrix(pbmc@assays$RNA@data[,batch1_index])
        g1=GROUP[batch1_index]
     
        while(j<=length(UBATCH)){
            batch2=UBATCH[j]
            batch2_index=which(BATCH==batch2)
            exp2=as.matrix(pbmc@assays$RNA@data[,batch2_index])
            g2=GROUP[batch2_index]
            VP_OUT=.getValidpair(exp1, g1, exp2, g2, CPU, method='kendall', print_step)
            this_vp=VP_OUT$vp 
            VP=cbind(VP,t(this_vp))
            j=j+1}
        i=i+1
        }

    VP=t(VP)
    DR=pbmc@reductions$pca@cell.embeddings


    

    MAP=rep('NA',length(GROUP))
    MAP[which(GROUP %in% VP[,1])]='V1'
    MAP[which(GROUP %in% VP[,2])]='V2'
    pbmc@meta.data$map=MAP
    

    OUT=.evaluateNew(DR, GROUP, VP)

    RESULT=list()
    RESULT$seurat=pbmc
    RESULT$vp=VP
    RESULT$vpcor=VP_OUT$cor
    RESULT$cor=OUT$cor
    RESULT$pv=OUT$pv
    RESULT$fdr=OUT$fdr
    return(RESULT)
}








