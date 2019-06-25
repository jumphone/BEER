source('BEER.R')

#Load Demo Data (Oligodendroglioma, GSE70630)
#Download: https://sourceforge.net/projects/beergithub/files/

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
rm(D12);rm(D34);rm(D56);rm(D1234);rm(D123456)




#DATA
#BATCH
MAXBATCH=''
GNUM=30
PCNUM=50
GN=2000
CPU=4
MTTAG="^MT-"
REGBATCH=FALSE
print_step=10
SEED=123
N=2

set.seed( SEED)
    RESULT=list()
    library(Seurat)
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('BEER start!')
    print(Sys.time())
    DATA=DATA
    BATCH=BATCH
    GNUM=GNUM
    PCNUM=PCNUM
    MTTAG=MTTAG
    MAXBATCH=MAXBATCH
    UBATCH=unique(BATCH)
    REGBATCH=REGBATCH
    GN=GN
    N=N
    print_step=print_step
    
    if(!MAXBATCH %in% UBATCH){
        MAXBATCH=names(which(table(BATCH)==max(table(BATCH))))
        }
    print('Max batch (MAXBATCH) is:')
    print(MAXBATCH)
    print('Group number (GNUM) is:')
    print(GNUM)
    print('Varible gene number (GN) is:')
    print(GN)
    
    VARG=c()
    for(this_batch in UBATCH){
        this_pbmc=CreateSeuratObject(counts = DATA[,which(BATCH==this_batch)], min.cells = 0, 
                                 min.features = 0, project = this_batch)
        this_pbmc <- NormalizeData(object = this_pbmc, normalization.method = "LogNormalize", 
                           scale.factor = 10000)
        this_pbmc <- FindVariableFeatures(object = this_pbmc, selection.method = "vst", nfeatures = GN)  
        this_varg=VariableFeatures(object = this_pbmc)
        VARG=c(VARG, this_varg)
        }
    VARG=unique(VARG)

pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL") 
    pbmc@meta.data$batch=BATCH
    VariableFeatures(object = pbmc)=VARG


    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = MTTAG)

    CPU=4
    if(REGBATCH==FALSE){
    pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc), vars.to.regress = c("nCount_RNA","percent.mt"), num.cores=CPU, do.par=TRUE)
    }else{
    pbmc <- ScaleData(object = pbmc, features = VariableFeatures(object = pbmc), vars.to.regress = c("nCount_RNA", "batch", "percent.mt"), num.cores=CPU, do.par=TRUE)
    }
    
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = VariableFeatures(object = pbmc), ndims.print=1,nfeatures.print=1)
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM,seed.use = SEED,n.components=N)




DR=pbmc@reductions$umap@cell.embeddings
    GROUP=rep('NA',length(BATCH))
    for(this_batch in UBATCH){
        this_index=which(BATCH==this_batch)
        this_one=DR[this_index,]
        #CNUM=max(c(5, round(length(this_index)/GNUM) ))
        this_group=.getGroup(this_one,this_batch,GNUM)
        GROUP[this_index]=this_group
    }


    pbmc@meta.data$group=GROUP






################
REF=.generate_ref(pbmc@assays$RNA@data, cbind(pbmc@meta.data$group,pbmc@meta.data$group),min_cell=1)

VREF=REF#REF[which(rownames(REF) %in% VARG),]
CVREF=cor(VREF,method='spearman')

SAME=3
VP=c()

I=1

library(igraph)
while(I<=TIME){
        
    p1=c()
    p2=c()
    score=c()
    i=1
    while(i<=nrow(CVREF)){
        this_p1=rownames(CVREF)[i]
        j=i+1
        while(j<=ncol(CVREF)){
            this_p2=colnames(CVREF)[j]  
            p1=c(p1,this_p1)
            p2=c(p2,this_p2)
            
            vp_index=which(VP[,1]==this_p1 & VP[,2]==this_p2) 
            if(length(vp_index) >0){
                this_score=999999 }else{
                this_score=1-CVREF[i,j]}         
            score=c(score,this_score)
            j=j+1}        
        i=i+1}


    NET = cbind(p1,p2) 
    g <- make_graph(t(NET),directed = FALSE)
    MST=mst(g, weights = score, algorithm = NULL)
    E_MST=as_edgelist(MST, names = TRUE)

    i=1
    while(i<=nrow(E_MST)){
        t1=unlist(strsplit(E_MST[i,1],'_'))[1]
        t2=unlist(strsplit(E_MST[i,2],'_'))[1]
        if(t1!=t2){VP=cbind(VP,E_MST[i,])}
            i=i+1}

    print(I)
    I=I+1
    }       
        
        
VP=t(VP)

#######################










DR=pbmc@reductions$pca@cell.embeddings  

    MAP=rep('NA',length(GROUP))
    MAP[which(GROUP %in% VP[,1])]='V1'
    MAP[which(GROUP %in% VP[,2])]='V2'
    pbmc@meta.data$map=MAP
    

    OUT=.evaluateProBEER(DR, GROUP, VP)


RESULT=list()
    RESULT$seurat=pbmc
    RESULT$vp=VP
    #RESULT$mst=MST
    RESULT$cor=OUT$cor
    RESULT$pv=OUT$pv
    RESULT$fdr=OUT$fdr
    RESULT$lcor=OUT$lcor
    RESULT$lpv=OUT$lpv
    RESULT$lc1=OUT$lc1
    RESULT$lc2=OUT$lc2
    RESULT$lfdr=OUT$lfdr
    
    PCUSE=which( (rank(RESULT$cor)>=length(RESULT$cor)/2 | RESULT$cor>0.7 )    & 
                (rank(RESULT$lcor) >=length(RESULT$cor)/2 | RESULT$lcor>0.7)   #&
                #p.adjust(RESULT$lc1,method='fdr') >0.05
               ) 
    
    RESULT$select=PCUSE



mybeer=RESULT
# Check selected PCs
PCUSE=mybeer$select
COL=rep('black',length(mybeer$cor))
COL[PCUSE]='red'
plot(mybeer$cor,mybeer$lcor,pch=16,col=COL,xlab='Rank Correlation',ylab='Linear Correlation',xlim=c(0,1),ylim=c(0,1))


pbmc <- mybeer$seurat
PCUSE <- mybeer$select
pbmc <- RunUMAP(object = pbmc, reduction.use='pca',dims = PCUSE, check_duplicates=FALSE)

DimPlot(pbmc, reduction.use='umap', group.by='batch', pt.size=0.1)

FeaturePlot(pbmc,features=c('AQP4','PDGFRA','NES','GFAP','TOP2A'))


