# R code to run ESTIMATE  to infer tumor purity for all samples in microarray.RData

require(estimate)
require(reshape2)
require(plyr)
require(dplyr)
load('server-side/RData/microarray.RData')

affy.probe.2.gene.name <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/affy.probe.2.gene.name.txt", stringsAsFactors=FALSE)
flag                   <- affy.probe.2.gene.name$HGNC.symbol == ''
affy.probe.2.gene.name <- affy.probe.2.gene.name[!flag,]

tmp                              <- ddply(affy.probe.2.gene.name,.(AFFY.HG.U133A.probe),nrow)
unambiguous.probe                <- tmp$AFFY.HG.U133A.probe[tmp$V1 == 1] %>% as.character #  keep probes mapping to only one gene
affy.probe.2.gene.name           <- affy.probe.2.gene.name[affy.probe.2.gene.name$AFFY.HG.U133A.probe %in% unambiguous.probe,]
rownames(affy.probe.2.gene.name) <- affy.probe.2.gene.name$AFFY.HG.U133A.probe

tmp                    <- ddply(affy.probe.2.gene.name,.(HGNC.symbol),nrow)
unambiguous.gene       <- tmp$HGNC.symbol[tmp$V1 == 1] %>% as.character # genes with only one probe mapped
expr.df                <- metastatic.breast.cancer.microarray.expr.matrix %>% as.data.frame
expr.df                <- expr.df[rownames(expr.df) %in% affy.probe.2.gene.name$AFFY.HG.U133A.probe,]
expr.df$HGNC.symbol    <- affy.probe.2.gene.name[rownames(expr.df),'HGNC.symbol']

merge.by.max <- function(x){
    x$HGNC.symbol <- NULL
    if(nrow(x) == 1){
        x
    }
    x$HGNC.symbol <- NULL
    m             <- as.matrix(x)
    value         <- apply(m,2,max)
    l             <- as.list(value)
    names(l)      <- colnames(x)
    as.data.frame(l)
}

merge.by.median <- function(x){
  x$HGNC.symbol <- NULL
  if(nrow(x) == 1){
    x
  }
  x$HGNC.symbol <- NULL
  m             <- as.matrix(x)
  value         <- apply(m,2,median)
  l             <- as.list(value)
  names(l)      <- colnames(x)
  as.data.frame(l)
}

merge.by.sum <- function(x){
  x$HGNC.symbol <- NULL
  if(nrow(x) == 1){
    x
  }
  x$HGNC.symbol <- NULL
  m             <- as.matrix(x)
  value         <- apply(m,2,sum)
  l             <- as.list(value)
  names(l)      <- colnames(x)
  as.data.frame(l)
}




unambiguous.only.df <- expr.df[expr.df$HGNC.symbol   %in% unambiguous.gene,] # expression profile of unambiguous genes
ambiguous.only.df   <- expr.df[!(expr.df$HGNC.symbol %in% unambiguous.gene),]
median.df           <- ddply(ambiguous.only.df,.(HGNC.symbol),merge.by.median)# expression profile of  all genes, probesets mapping to the same gene were combined by median 
median.df           <- rbind(median.df,unambiguous.only.df)
max.df              <- ddply(ambiguous.only.df,.(HGNC.symbol),merge.by.max)# expression profile of unambiguous genes, probesets mapping to the same gene were combined by max
max.df              <- rbind(max.df,unambiguous.only.df)
sum.df              <- ddply(ambiguous.only.df,.(HGNC.symbol),merge.by.sum)# expression profile of unambiguous genes, probesets mapping to the same gene were combined by sum
sum.df              <- rbind(sum.df,unambiguous.only.df)




########## run estimate###############
data(common_genes)

setwd('~/Project/Cancer2CellLine/client-side/')
median.df             <- median.df[median.df$HGNC.symbol %in% common_genes$GeneSymbol,]
rownames(median.df)   <- median.df$HGNC.symbol
median.df$HGNC.symbol <- NULL
outputGCT(input.f = median.df,output.f = 'output/run.estimate.R.output/median.df.gct')
estimateScore(input.ds = 'output/run.estimate.R.output/median.df.gct',output.ds = 'output/run.estimate.R.output/median.df.score.gct',platform = 'affymetrix')

max.df             <- max.df[max.df$HGNC.symbol %in% common_genes$GeneSymbol,]
rownames(max.df)   <- max.df$HGNC.symbol
max.df$HGNC.symbol <- NULL
outputGCT(input.f = max.df,output.f = 'output/run.estimate.R.output/max.df.gct')
estimateScore(input.ds = 'output/run.estimate.R.output/max.df.gct',output.ds = 'output/run.estimate.R.output/max.df.score.gct',platform = 'affymetrix')

sum.df             <- sum.df[sum.df$HGNC.symbol %in% common_genes$GeneSymbol,]
rownames(sum.df)   <- sum.df$HGNC.symbol
sum.df$HGNC.symbol <- NULL
outputGCT(input.f = sum.df,output.f = 'output/run.estimate.R.output/sum.df.gct')
estimateScore(input.ds = 'output/run.estimate.R.output/sum.df.gct',output.ds = 'output/run.estimate.R.output/sum.df.score.gct',platform = 'affymetrix')

unambiguous.only.df             <- unambiguous.only.df[unambiguous.only.df$HGNC.symbol %in% common_genes$GeneSymbol,]
rownames(unambiguous.only.df)   <- unambiguous.only.df$HGNC.symbol
unambiguous.only.df$HGNC.symbol <- NULL
outputGCT(input.f = unambiguous.only.df,output.f = 'output/run.estimate.R.output/unambiguous.only.df.gct')
estimateScore(input.ds = 'output/run.estimate.R.output/unambiguous.only.df.gct',output.ds = 'output/run.estimate.R.output/unambiguous.only.df.score.gct',platform = 'affymetrix')

setwd('~/Project/Cancer2CellLine/')

require(data.table)
max.df.score              <- fread('./client-side/output/run.estimate.R.output//max.df.score.gct',header=TRUE,skip=2)  %>% as.data.frame
median.df.score           <- fread('./client-side/output/run.estimate.R.output//median.df.score.gct',header=TRUE,skip=2) %>% as.data.frame
sum.df.score              <- fread('./client-side/output/run.estimate.R.output//sum.df.score.gct',header=TRUE,skip=2) %>% as.data.frame
unambiguous.only.df.score <- fread('./client-side/output/run.estimate.R.output//unambiguous.only.df.score.gct',header=TRUE,skip=2) %>% as.data.frame

m <- cbind(max.df.score[4,3:119]  %>% unlist ,median.df.score[4,3:119] %>% unlist,sum.df.score[4,3:119] %>% unlist,unambiguous.only.df.score[4,3:119] %>% unlist)
pairs(m)


load('server-side/RData/microarray.RData')
biopsy.site       <-  metastatic.breast.cancer.microarray.meta.df$biopsy.site
liver.samples     <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'LIVER']  %>% as.character
lung.samples      <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'LUNG']   %>% as.character
bone.samples      <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'BONE']   %>% as.character
brain.samples     <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'BRAIN']  %>% as.character
lymph.samples     <-  metastatic.breast.cancer.microarray.meta.df$gsm.id[biopsy.site == 'LYMPH_NODE']  %>% as.character


biopsy.site       <-  metastatic.breast.cancer.microarray.meta.df[match(rownames(m) %>% as.character,metastatic.breast.cancer.microarray.meta.df$gsm.id),'biopsy.site']
df                <-  data.frame(biopsy.site = biopsy.site,estimate.score = m[,2]) # OK, for genes with multiple probes, use the median value as its expression
df                <-  df[df$biopsy.site %in% c('LIVER','LUNG','LYMPH_NODE','BONE','BRAIN'),]
microarray.data.tumor.purity.df <- df

save(file = 'client-side/output/run.estimate.R.output/microarray.data.tumor.purity.RData',list=c('microarray.data.tumor.purity.df'))
