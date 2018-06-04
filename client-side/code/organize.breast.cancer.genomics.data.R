### R code to organize somatic mutation and CNV data for MET500,CCLE and TCGA ###

require(plyr)
require(dplyr)
require(reshape2)
require(CNTools)
require(cgdsr)
source('client-side/code/util.R')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')


######     Let us first solve the gene symbol issues. It is really nasty! ############
########## Finally, map all CCLE captured gene symbols to the newest HGNC version   ################
hgnc.complete.set               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/hgnc_complete_set.txt", stringsAsFactors=FALSE)
hgnc.complete.set               <- subset( hgnc.complete.set,select=c('hgnc_id','symbol'))
colnames(hgnc.complete.set)     <- c('HGNC.ID','Symbol')
df                              <- read.delim("~/Project/Cancer2CellLine/client-side/Data/somatic.mutation/CCLE.captured.gene.txt", stringsAsFactors=FALSE)
df                              <- df[df$HGNC.Symbol != '',]
idx                             <- match(x = df$HGNC.ID,table = hgnc.complete.set$HGNC.ID)
mapping.df                      <- data.frame(previous.hgnc.symbol=df$HGNC.Symbol,new.hgnc.symbol=hgnc.complete.set$Symbol[idx]) %>% unique # There are duplicates in the CCLE.captured.gene.txt: two genes correspond to HGNC symbol MECOM
rownames(mapping.df)            <- mapping.df$previous.hgnc.symbol %>% as.character
mapping.df$previous.hgnc.symbol <- as.character(mapping.df$previous.hgnc.symbol)
mapping.df$new.hgnc.symbol      <- as.character(mapping.df$new.hgnc.symbol) # The mapping.df is used to map all HGNC symbols to newest version


hg19.RefSeq.gene.coordinates <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/hg19.RefSeq.gene.coordinates.txt", stringsAsFactors=FALSE)
flag                         <- grepl(x=hg19.RefSeq.gene.coordinates$chrom,pattern = '_')
hg19.RefSeq.gene.coordinates <- hg19.RefSeq.gene.coordinates[!flag,] #Let us remove the genes which are on the un-assembled contigs!
tmp                          <- ddply(hg19.RefSeq.gene.coordinates,.(name2),function(x) data.frame(chrom=x$chrom[1], start=min(x$txStart),end=max(x$txEnd),geneid=x$name2[1],genename=x$name2[1]))
tmp$name2                    <- NULL
hg19.gene.info               <- tmp
hg19.gene.info$chrom         <- gsub(x=hg19.gene.info$chrom,pattern='chr',replacement = '')
hg19.gene.info$chrom         <- as.character(hg19.gene.info$chrom)
hg19.gene.info$geneid        <- as.character(hg19.gene.info$geneid)
hg19.gene.info$genename      <- as.character(hg19.gene.info$genename)
flag                         <- which(hg19.gene.info$genename %in% mapping.df$previous.hgnc.symbol)
for(i in flag){
    pre.symbol <- hg19.gene.info$genename[i]  %>% as.character
    new.symbol <- mapping.df[pre.symbol,'new.hgnc.symbol']
    hg19.gene.info[i,'genename']        <- new.symbol
    hg19.gene.info[i,'geneid']          <- new.symbol
} #In hg19 RefSeq release, some CCLE captured genes are still with expired HGNC symbol, just replace them

CCLE.captured.gene.hg19.info       <- hg19.gene.info[hg19.gene.info$genename %in% mapping.df$new.hgnc.symbol,] # Well, I trust RefSeq. If some gene symbols here were missed after this filtering, then remove them. Actually, only 8 were missed 
CCLE.captured.gene.hg19.info       <- CCLE.captured.gene.hg19.info[CCLE.captured.gene.hg19.info$chrom != 'Y',] # Well, we are studying breast cancer, remove genes on chromosome Y. 
CCLE.captured.gene                 <- CCLE.captured.gene.hg19.info$genename



################ Organize CCLE gene mutation data ################
####oject to be saved: 
## CCLE.breast.cancer.gene.mutation.profile 
## CCLE.breast.cancer.gene.mutation.freq########

CCLE.maf.data               <- read.delim("~/Project/Cancer2CellLine/client-side/Data/somatic.mutation/CCLE_hybrid_capture1650_hg19_NoCommonSNPs_NoNeutralVariants_CDS_2012.05.07.maf",stringsAsFactors=FALSE)
flag                        <- CCLE.maf.data$Hugo_Symbol %in% mapping.df$previous.hgnc.symbol # In the CCLE maf files,there are some genes outside the capture-seq probe regions. Currently I just remove them,only consider CCLE captured genes
CCLE.maf.data               <- CCLE.maf.data[flag,]
idx                         <- match(CCLE.maf.data$Hugo_Symbol,mapping.df$previous.hgnc.symbol)
CCLE.maf.data$Hugo_Symbol   <- mapping.df$new.hgnc.symbol[idx] # map previous HGNC symbols to newest version 
flag                        <- CCLE.maf.data$Hugo_Symbol %in% CCLE.captured.gene
CCLE.maf.data               <- CCLE.maf.data[flag,]
flag                        <- CCLE.maf.data$Tumor_Sample_Barcode %in% CCLE.genotyped.and.rna.seqqed.cell.line
CCLE.maf.data               <- CCLE.maf.data[flag,]
flag                        <- CCLE.maf.data$Tumor_Sample_Barcode %in% CCLE.breast.cancer.cell.line
CCLE.breast.cancer.maf.data <- CCLE.maf.data[flag,]


CCLE.gene.mutation.profile           <- matrix(data=0,nrow=CCLE.captured.gene %>% length,ncol=CCLE.genotyped.and.rna.seqqed.cell.line %>% length)
rownames(CCLE.gene.mutation.profile) <- CCLE.captured.gene
colnames(CCLE.gene.mutation.profile) <- CCLE.genotyped.and.rna.seqqed.cell.line
for (i in 1:nrow(CCLE.maf.data)) {
    gene      <- CCLE.maf.data[i,'Hugo_Symbol']          %>% as.character
    cell.line <- CCLE.maf.data[i,'Tumor_Sample_Barcode'] %>% as.character
    CCLE.gene.mutation.profile[gene,cell.line] <- 1
}

CCLE.breast.cancer.gene.mutation.profile <- CCLE.gene.mutation.profile[,CCLE.breast.cancer.cell.line]
CCLE.breast.cancer.gene.mutation.freq    <- apply(CCLE.breast.cancer.gene.mutation.profile,1,function(x) sum(x)/ncol(CCLE.breast.cancer.gene.mutation.profile)) %>% sort(decreasing = TRUE)



############ organize MET500 gene mutation data ##################
#objects to be saved: 
## MET500.breast.cancer.gene.mutation.profile 
## MET500.breast.cancer.gene.mutation.freq
## MET500.breast.cancer.top.mutated.gene

load('server-side/RData/MET500.RData')
get.id <- function(x){
    l <- strsplit(x = x,split='\\.')  %>% unlist 
    l[1]
}
flag                           <- grepl(x=MET500.sample.meta$cancer.type,pattern='Breast Invasive Ductal Carcinoma') 
breast.cancer.sample.MET500.id <- MET500.sample.meta$MET500.id[flag] %>% as.character %>% unique

MET500.maf.data                <- read.csv("~/Project/Cancer2CellLine/client-side/Data/somatic.mutation/somatic_v4.csv",stringsAsFactors=FALSE)
MET500.maf.data$MET500.id      <- sapply(MET500.maf.data$Pipeline_ID,get.id)
MET500.maf.data                <- MET500.maf.data[,c('MET500.id','Gene','Effect')] 

# Well, let us only care about the gene symbols which are presented in CCLE. There might be also expired gene symboles here 
flag            <- MET500.maf.data$Gene %in% mapping.df$previous.hgnc.symbol | MET500.maf.data$Gene %in% mapping.df$new.hgnc.symbol
MET500.maf.data <- MET500.maf.data[flag,]
flag            <- which(MET500.maf.data$Gene %in% mapping.df$previous.hgnc.symbol)
for(i in flag){
    pre.symbol <- MET500.maf.data$Gene[i]  %>% as.character
    new.symbol <- mapping.df[pre.symbol,'new.hgnc.symbol']
    MET500.maf.data[i,'Gene']        <- new.symbol
} # map previous symbols to newest symbol

MET500.breast.cancer.maf.data                        <- MET500.maf.data[MET500.maf.data$MET500.id %in% breast.cancer.sample.MET500.id,]
MET500.breast.cancer.gene.mutation.profile           <- matrix(data=0,nrow=CCLE.captured.gene %>% length,ncol = breast.cancer.sample.MET500.id %>% length)
rownames(MET500.breast.cancer.gene.mutation.profile) <- CCLE.captured.gene
colnames(MET500.breast.cancer.gene.mutation.profile) <- breast.cancer.sample.MET500.id
for (i in 1:nrow(MET500.breast.cancer.maf.data)) {
    gene      <- MET500.breast.cancer.maf.data[i,'Gene']      %>% as.character
    sample.id <- MET500.breast.cancer.maf.data[i,'MET500.id'] %>% as.character
    MET500.breast.cancer.gene.mutation.profile[gene,sample.id] <- 1
}
MET500.breast.cancer.gene.mutation.freq    <- apply(MET500.breast.cancer.gene.mutation.profile,1,function(x) sum(x)/ncol(MET500.breast.cancer.gene.mutation.profile)) %>% sort(decreasing = TRUE)
MET500.breast.cancer.top.mutated.gene      <- names(MET500.breast.cancer.gene.mutation.freq)[MET500.breast.cancer.gene.mutation.freq >= 0.05]




################ Organize TCGA mutation data, I only download the mutation data for CCLE captured genes #######
#Objects to be saved:
##TCGA.breast.cancer.gene.mutation.profile
##TCGA.breast.cancer.gene.mutation.freq

#TCGA.breast.cancer.top.mutated.gene <- c('TP53','PIK3CA','TTN','MAP3K1','KMT2C') # hotspot genes showing mutation freq >= 0.05 in Breast Invasive Ductal Carcinoma 
#GATA3,MUC16,MUC12 is not captured in CCLE. Access http://www.cbioportal.org/study?id=brca_tcga#summary, make sure that in the cancer type detailed tab Breast Invasive Ductal Carcinoma  is selected
#pool.gene             <- c(TCGA.breast.cancer.top.mutated.gene,MET500.breast.cancer.top.mutated.gene) %>% unique

study.id   <- 'brca_tcga'
mycgds     <-  CGDS("http://www.cbioportal.org/public-portal/")
case.list  <- 'brca_tcga_sequenced'
profile.id <- 'brca_tcga_mutations'
clinical.data                           <- getClinicalData(mycgds,case.list)
Breast.Invasive.Ductal.Carcinoma.sample <- rownames(clinical.data)[clinical.data$CANCER_TYPE_DETAILED == 'Breast Invasive Ductal Carcinoma']

tmp3 <- foreach(i = seq(from = 1,to=length(CCLE.captured.gene)-10,by=10)) %do% {
    low <- i
    up  <- i + 9
    if(up > length(CCLE.captured.gene) ) {
        up <- length(CCLE.captured.gene)
    }
    pool.gene <- CCLE.captured.gene[low:up]
    tmp      <- getProfileData(mycgds,pool.gene,profile.id,case.list)
    tmp      <- tmp[rownames(tmp) %in% Breast.Invasive.Ductal.Carcinoma.sample,]
    TCGA.breast.cancer.gene.mutation.profile <- as.matrix(tmp) 
    TCGA.breast.cancer.gene.mutation.profile[is.na(TCGA.breast.cancer.gene.mutation.profile)]    <- '0'
    TCGA.breast.cancer.gene.mutation.profile[TCGA.breast.cancer.gene.mutation.profile == 'NaN' ] <- '0'
    TCGA.breast.cancer.gene.mutation.profile[TCGA.breast.cancer.gene.mutation.profile != '0']    <- '1'
    tmp2 <- matrix(data=as.integer(TCGA.breast.cancer.gene.mutation.profile), 
                   ncol=ncol(TCGA.breast.cancer.gene.mutation.profile), 
                   nrow=nrow(TCGA.breast.cancer.gene.mutation.profile)
                   )
    rownames(tmp2) <- rownames(TCGA.breast.cancer.gene.mutation.profile)
    colnames(tmp2) <- colnames(TCGA.breast.cancer.gene.mutation.profile)
    tmp2
}

row.name <- rownames(tmp3[[1]])
tmp4 <- foreach(l = tmp3,.combine='cbind') %do% {
    l[row.name,]  
}


TCGA.breast.cancer.gene.mutation.profile <- tmp4 %>% t
flag                                     <- rownames(TCGA.breast.cancer.gene.mutation.profile) == 'NKX2.1'
rownames(TCGA.breast.cancer.gene.mutation.profile)[flag] <- 'NKX2-1' #hmm, TCGA returns NKX2-1 for gene NKX2.1, replace back
TCGA.breast.cancer.gene.mutation.freq    <- apply(TCGA.breast.cancer.gene.mutation.profile,1,function(x) sum(x)/ncol(TCGA.breast.cancer.gene.mutation.profile)) %>% sort(decreasing = TRUE)
TCGA.breast.cancer.top.mutated.gene      <- names(TCGA.breast.cancer.gene.mutation.freq)[TCGA.breast.cancer.gene.mutation.freq >= 0.05]

######################### Save the data into RData file format #############
save(file='client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.mutation.profile.RData',
     list=c('CCLE.breast.cancer.gene.mutation.profile',
            'CCLE.breast.cancer.gene.mutation.freq',
            'CCLE.gene.mutation.profile',
            'MET500.breast.cancer.gene.mutation.profile',
            'MET500.breast.cancer.gene.mutation.freq',
            'MET500.breast.cancer.top.mutated.gene',
            'TCGA.breast.cancer.gene.mutation.profile',
            'TCGA.breast.cancer.gene.mutation.freq',
            'TCGA.breast.cancer.top.mutated.gene'
            )
    )



########## Organize CCLE cnv data, with CNTools package #########
##objects to be saved:
#CCLE.breast.cancer.gene.cnv
#CCLE.breast.cancer.gene.cnv.median
#CCLE.gene.cnv

CCLE.cnv.data               <- read.delim("~/Project/Cancer2CellLine/client-side/Data/CNV/CCLE_copynumber_2013-12-03.seg.txt", stringsAsFactors=FALSE)
colnames(CCLE.cnv.data)     <- c('ID','chrom','loc.start','loc.end','num.mark','seg.mean')
CCLE.seg                    <- CNSeg(CCLE.cnv.data)
CCLE.gene.cnv               <- getRS(object = CCLE.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = CCLE.captured.gene.hg19.info,what = "max")@rs
CCLE.gene.cnv               <- CCLE.gene.cnv[,5:ncol(CCLE.gene.cnv)]
rownames(CCLE.gene.cnv)     <- CCLE.gene.cnv$genename
CCLE.gene.cnv$genename      <- NULL
CCLE.gene.cnv               <- as.matrix(CCLE.gene.cnv)
flag                        <- grepl(x = colnames(CCLE.gene.cnv),pattern = 'BREAST')
CCLE.breast.cancer.gene.cnv <- CCLE.gene.cnv[,flag]
CCLE.breast.cancer.gene.cnv.median    <- apply(CCLE.breast.cancer.gene.cnv,1,median) 

CCLE.breast.cancer.gene.cnv.median.metastatic     <- apply(CCLE.breast.cancer.gene.cnv[,CCLE.breast.cancer.metastatic.cell.line],1,median)
CCLE.breast.cancer.gene.cnv.median.non.metastatic <- apply(CCLE.breast.cancer.gene.cnv[,CCLE.breast.cancer.non.metastatic.cell.line],1,median)




####### Organize MET500 gene cnv data ###########
##objects to be saved:
#MET500.breast.cancer.gene.cnv
#MET500.breast.cancer.gene.cnv.median
MET500.cnv.data               <- read.csv("~/Project/Cancer2CellLine/client-side/Data/CNV/cnv_v4.csv", stringsAsFactors=FALSE)
MET500.cnv.data               <- MET500.cnv.data[,c('Pipeline_ID','Chr','Start','End','Log2_Coverage_Ratio')]
colnames(MET500.cnv.data)     <- c('ID','chrom','loc.start','loc.end','seg.mean')
MET500.cnv.data$ID            <- sapply(MET500.cnv.data$ID,get.id)
MET500.breast.cancer.cnv.data <- MET500.cnv.data[MET500.cnv.data$ID %in% breast.cancer.sample.MET500.id,]
MET500.breast.cancer.seg      <- CNSeg(MET500.breast.cancer.cnv.data)
MET500.breast.cancer.gene.cnv <- getRS(object = MET500.breast.cancer.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = CCLE.captured.gene.hg19.info,what = "max")@rs
MET500.breast.cancer.gene.cnv <- MET500.breast.cancer.gene.cnv[,5:ncol(MET500.breast.cancer.gene.cnv)]
rownames(MET500.breast.cancer.gene.cnv) <- MET500.breast.cancer.gene.cnv$genename
MET500.breast.cancer.gene.cnv$genename  <- NULL
MET500.breast.cancer.gene.cnv           <- as.matrix(MET500.breast.cancer.gene.cnv)
MET500.breast.cancer.gene.cnv.median    <- apply(MET500.breast.cancer.gene.cnv,1,median) 


###### Organize TCGA cnv data ############
##Objects to be saved:
#TCGA.breast.cancer.gene.cnv
#TCGA.breast.cancer.gene.cnv.median

TCGA.cnv.data <- read.delim("~/Project/Cancer2CellLine/client-side/Data/CNV/BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt", stringsAsFactors=FALSE)
convert.id.style <- function(x){
   l    <- strsplit(x,split = '-') %>% unlist  
   l[4] <- gsub(x=l[4],pattern='[:A-Z:]',replacement='')
   paste(l[1:4],sep = '',collapse = '.')
}
TCGA.cnv.data$Sample        <- sapply(TCGA.cnv.data$Sample,convert.id.style)
flag                        <- TCGA.cnv.data$Sample %in% Breast.Invasive.Ductal.Carcinoma.sample
TCGA.breast.cancer.cnv.data <- TCGA.cnv.data[flag,]
colnames(TCGA.breast.cancer.cnv.data) <- c('ID','chrom','loc.start','loc.end','num.marker','seg.mean')
TCGA.breast.cancer.seg      <- CNSeg(TCGA.breast.cancer.cnv.data)
TCGA.breast.cancer.gene.cnv <- getRS(object = TCGA.breast.cancer.seg,by = 'gene',imput = FALSE,XY=TRUE,geneMap = CCLE.captured.gene.hg19.info,what = "max")@rs
TCGA.breast.cancer.gene.cnv <- TCGA.breast.cancer.gene.cnv[,5:ncol(TCGA.breast.cancer.gene.cnv)]
rownames(TCGA.breast.cancer.gene.cnv) <- TCGA.breast.cancer.gene.cnv$genename
TCGA.breast.cancer.gene.cnv$genename  <- NULL
TCGA.breast.cancer.gene.cnv           <- as.matrix(TCGA.breast.cancer.gene.cnv)
TCGA.breast.cancer.gene.cnv.median    <- apply(TCGA.breast.cancer.gene.cnv,1,median) 

######################### Save the data into RData file format #############
save(file='client-side/output/organize.breast.cancer.genomics.data.R.output/breast.cancer.gene.cnv.RData',
     list=c('CCLE.breast.cancer.gene.cnv',
            'CCLE.breast.cancer.gene.cnv.median','CCLE.breast.cancer.gene.cnv.median.metastatic','CCLE.breast.cancer.gene.cnv.median.non.metastatic',
            'MET500.breast.cancer.gene.cnv',
            'MET500.breast.cancer.gene.cnv.median',
            'TCGA.breast.cancer.gene.cnv',
            'TCGA.breast.cancer.gene.cnv.median',
            'CCLE.gene.cnv'
     )
)

