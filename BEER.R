# Batch EffEct Remover for single-cell data (BEER)
# Author: Feng Zhang
# Date: Mar. 7, 2019

#source('https://raw.githubusercontent.com/jumphone/BEER/master/BEER.R')

#library(Seurat)
#library(pcaPP)


############################################################################################
############################################################################################

.get_cor  <- function(exp_sc_mat, exp_ref_mat, method='kendall',CPU=4, print_step=10, gene_check=FALSE){
    #method = "pearson", "kendall", "spearman"
    ##################
    print('Gene number of exp_sc_mat1:')
    print(nrow(exp_sc_mat))
    print('Gene number of exp_sc_mat2:')
    print(nrow(exp_ref_mat))
    #################
    library(parallel)
    #Step 1. get overlapped genes
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    ###############
    print('Number of overlapped genes:')
    print(nrow(exp_sc_mat))
    if(gene_check==TRUE){
    print('Press RETURN to continue:')
    scan();}
    ###################
    #Step 2. calculate prob
    SINGLE <- function(i){
        library('pcaPP')
        exp_sc = as.array(exp_sc_mat[,i])
        
        cor_list=rep(0,length(colname_ref))
        j=1
        while(j<=length(colname_ref)){
            exp_ref = as.array(exp_ref_mat[,j])
            
            if(method=='kendall'){this_cor=cor.fk(exp_sc,exp_ref)}
            else{
            this_cor=cor(exp_sc,exp_ref, method=method)}
            
            cor_list[j]=this_cor
            j=j+1}
        ################################
        if(i%%print_step==1){print(i)}
        return(cor_list)
        }
    #######################################
    cl= makeCluster(CPU,outfile='')
    RUN = parLapply(cl=cl,1:length(exp_sc_mat[1,]), SINGLE)
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    COR = c()
    for(cor_list in RUN){
        COR=cbind(COR, cor_list)}
    #######################################
    rownames(COR)=colname_ref
    colnames(COR)=colname_sc
    return(COR)
    }


.get_tag_max <- function(COR){
    RN=rownames(COR)
    CN=colnames(COR)
    TAG=cbind(CN,rep('NA',length(CN)))
    i=1
    while(i<=length(CN)){
        this_rn_index=which(COR[,i] == max(COR[,i]))[1]
        TAG[i,2]=RN[this_rn_index]
        i=i+1
        }
    colnames(TAG)=c('cell_id','tag')
    return(TAG)
    }

.simple_combine <- function(exp_sc_mat1, exp_sc_mat2){    
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }
    


.generate_ref <- function(exp_sc_mat, TAG, min_cell=1, refnames=FALSE){
    NewRef=c()
    TAG[,2]=as.character(TAG[,2])
    if(refnames==FALSE){
        refnames=names(table(TAG[,2]))}
        else{refnames=refnames}
    outnames=c()
    for(one in refnames){
        this_col=which(TAG[,2]==one)
        if(length(this_col)>= min_cell){
            outnames=c(outnames,one)
            if(length(this_col) >1){
                this_new_ref=apply(exp_sc_mat[,this_col],1,sum)
                }
                else{this_new_ref = exp_sc_mat[,this_col]}
            NewRef=cbind(NewRef,this_new_ref)
            }
        }
    rownames(NewRef)=rownames(exp_sc_mat)
    colnames(NewRef)=outnames
    if(length(NewRef[1,])==1){
        NewRef=cbind(NewRef[,1], NewRef[,1])
        rownames(NewRef)=rownames(exp_sc_mat)
        colnames(NewRef)=c(outnames,outnames)
        }
    return(NewRef)
    }



############################################################################################
############################################################################################



.data2one <- function(DATA, GENE, CPU=4, PCNUM=100){
    PCUSE=1:PCNUM
    print('Start')
    library(Seurat)
    print('Step1.Create Seurat Object...')
    DATA = CreateSeuratObject(raw.data = DATA, min.cells = 0, min.genes = 0, project = "DATA") 
    print('Step2.Normalize Data...')
    DATA <- NormalizeData(object = DATA, normalization.method = "LogNormalize", scale.factor = 10000)
    print('Step3.Scale Data...')
    DATA <- ScaleData(object = DATA, genes.use =GENE, vars.to.regress = c("nUMI"), num.cores=CPU, do.par=TRUE)
    print('Step4.PCA...')
    DATA <- RunPCA(object = DATA, pcs.compute=PCNUM, pc.genes =GENE, do.print = FALSE)
    print('Step5.One-dimention...')
    DATA <- RunTSNE(object = DATA, dims.use = PCUSE, do.fast=TRUE,dim.embed = 1)
    DR=DATA@dr$tsne@cell.embeddings
    print('Finished!!!')
    return(DR)
    }

.getGroup <- function(X,TAG,CNUM=100){
    DR=X
    RANK=rank(DR,ties.method='random')
    CUTOFF=CNUM 
    GROUP=rep('NA',length(RANK))
    i=1
    j=1
    while(i<=length(RANK)){
        GROUP[which(RANK==i)]=paste0(TAG,'_',as.character(j))
        #if(i%%CUTOFF==1){j=j+1;print(j)}
        if(i%%CUTOFF==1){j=j+1}
        i=i+1}
    print('Group Number:')
    print(j-1)
    return(GROUP)
}


.getValidpair <- function(DATA1, GROUP1, DATA2, GROUP2, CPU=4, method='kendall', print_step=10){
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
    print('Start')
    print('Step1.Generate Reference...')
    REF1=.generate_ref(DATA1, cbind(GROUP1, GROUP1), min_cell=1) 
    REF2=.generate_ref(DATA2, cbind(GROUP2, GROUP2), min_cell=1) 
    print('Step2.Calculate Correlation Coefficient...')
    out = .get_cor( REF1, REF2, method=method,CPU=CPU, print_step=print_step)
    print('Step3.Analyze Result...')
    tag1=.get_tag_max(out)
    tag2=.get_tag_max(t(out))
    V=c()
    i=1
    while(i<=nrow(tag1)){
        t1=tag1[i,1]
        t2=tag1[i,2]
        if(tag2[which(tag2[,1]==t2),2]==t1){V=c(V,i)}           
        i=i+1}
    VP=tag1[V,]
    C=c()
    t=1
    while(t<=nrow(VP)){
        this_c=out[which(rownames(out)==VP[t,2]),which(colnames(out)==VP[t,1])]
        C=c(C,this_c)
        t=t+1}
    #if(do.plot==TRUE){plot(C)}
    #VP=VP[which(C>=CUTOFF),]  
    print('Finished!!!')
    OUT=list()
    OUT$vp=VP
    OUT$cor=C
    return(OUT)
    }



.evaluateBatcheffect <- function(DR, B1index, B2index, GROUP, VP){
    
    #library(dtw)
    OUT=list()
    
    VALID_PAIR=VP
    
    ALL_COR=c()   
    ALL_PV=c() 
    
    index1=B1index
    index2=B2index 
    
    print('Start')
    THIS_DR=1
    while(THIS_DR<=ncol(DR)){
        THIS_PC = DR[,THIS_DR]
        
        all_lst1=DR[index1,THIS_DR]
        all_lst2=DR[index2,THIS_DR] 
        
        maplst1=c()
        maplst2=c()

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




BEER <- function(D1, D2, CNUM=10, PCNUM=50, VPCOR=0, CPU=4, print_step=10){
    RESULT=list()
    library(Seurat)
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
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


# BEER with Multiple Samples
MBEER <- function(DATA, BATCH, MAXBATCH=NULL, CNUM=10, PCNUM=50, CPU=4, print_step=10){
  
    RESULT=list()
    library(Seurat)
    #source('https://raw.githubusercontent.com/jumphone/scRef/master/scRef.R')
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
    if(!MAXBATCH %in% rownames(TABLE)){
        MAXBATCH=rownames(TABLE)[which(TABLE==max(TABLE))[1]]}
    print(MAXBATCH)
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
    print('Analyze Pair:')
    print(i)
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
          
    i=i+1}
    
    colnames(COR)=PAIR[,2]
    colnames(PV)=PAIR[,2]
    colnames(FDR)=PAIR[,2]
       
    print('############################################################################')
    print('MainStep3. Output')
    print('############################################################################')
    
    ########################## 
    RESULT$seurat=pbmc
    RESULT$COR=COR
    RESULT$PV=PV
    RESULT$FDR=FDR
    RESULT$MAXBATCH=MAXBATCH
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
