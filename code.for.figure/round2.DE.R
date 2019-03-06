### R code to perform differential expression analysis between MET500 breast cancer samples , CCLE cell lines and organoids ###
### RUVg inferred hidden factors are plugged into model fitting

load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
load('client-side/output/pick.out.cell.line.by.subtype.R.output/pick.out.cell.line.by.subtype.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/organize.breast.cancer.organoid.data.R.output/organoid.RData')
require(DESeq2)
require(dplyr)
require(RUVSeq)

MET500.breast.cancer.polyA.non.Basal.sample <- c(MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.Her2.sample)
MET500.breast.cancer.polyA.Basal.sample     <- c(MET500.breast.cancer.polyA.Basal.sample)
CCLE.non.Basal.cell.line                    <- c(CCLE.breast.cancer.Her2.cell.line,CCLE.breast.cancer.LumA.cell.line,CCLE.breast.cancer.LumB.cell.line)
non.Basal.organoid                          <- c(LumA.organoid,LumB.organoid,Her2.organoid)

common.gene         <- intersect(rownames(CCLE.log2.rpkm.matrix),rownames(MET500.log2.fpkm.matrix))
common.gene         <- intersect(common.gene,rownames(organoid.log2.rpkm.matrix))
tmp                 <- cbind(CCLE.log2.read.count.matrix[common.gene,CCLE.breast.cancer.cell.line],organoid.log2.read.count.matrix[common.gene,])
tmp                 <- cbind(tmp,MET500.log2.read.count.matrix[common.gene,MET500.breast.cancer.polyA.sample])
flag                <- apply(tmp,1,function(x) sum(x>=1)) # filter lowly expressed genes
common.gene         <- common.gene[flag >=10]

biomart.protein.coding.gene <- read.table("~/Project/Cancer2CellLine/client-side/meta.data/biomart.protein.coding.gene.txt", quote="\"", stringsAsFactors=FALSE)$V1 %>% as.character
common.gene <- intersect(common.gene,biomart.protein.coding.gene)


############ DE analysis for non-Basal-cell-lines vs non-Basal-MET500 samples#########

expr.matrix  <- cbind(CCLE.log2.read.count.matrix[common.gene,CCLE.non.Basal.cell.line],MET500.log2.read.count.matrix[common.gene,MET500.breast.cancer.polyA.non.Basal.sample])
df           <- data.frame(condition=c(rep(x='cell.line',times=length(CCLE.non.Basal.cell.line)), rep(x='MET500',times=length(MET500.breast.cancer.polyA.non.Basal.sample))))
df$condition <- factor(df$condition,levels = c('cell.line','MET500'))
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition)

dds <- DESeq(dds )
res <- results(dds)
res <- res[complete.cases(res),]
res <- res[order(res$pvalue),]
de.gene   <- rownames(res)[res$padj < 0.01 ]



RUVg.rs <- RUVg(x = round(2^expr.matrix-1),cIdx = ifelse(rownames(expr.matrix) %in% de.gene,FALSE,TRUE),k=1)
W <- RUVg.rs$W
rownames(W) <- colnames(expr.matrix)
df$W <- W[,1]
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition + W)

dds <- DESeq(dds)
res <- results(dds,contrast=c('condition','MET500','cell.line'))
res <- res[complete.cases(res),]
res <- res[order(res$pvalue),]
non.Basal.MET500.vs.CCLE.down.de.gene <- (res)[res$padj < 0.001 & res$log2FoldChange < -1,]
non.Basal.MET500.vs.CCLE.up.de.gene   <- (res)[res$padj < 0.001 & res$log2FoldChange  > 1,]




############ DE analysis for Basal-cell-lines vs Basal-MET500 samples#########
CCLE.Basal.cell.line <- CCLE.breast.cancer.Basal.cell.line

expr.matrix  <- cbind(CCLE.log2.read.count.matrix[common.gene,CCLE.Basal.cell.line],MET500.log2.read.count.matrix[common.gene,MET500.breast.cancer.polyA.Basal.sample])
df           <- data.frame(condition=c(rep(x='cell.line',times=length(CCLE.Basal.cell.line)), rep(x='MET500',times=length(MET500.breast.cancer.polyA.Basal.sample))))
df$condition <- factor(df$condition,levels = c('cell.line','MET500'))
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res <- res[complete.cases(res),]
res <- res[order(res$pvalue),]
de.gene   <- rownames(res)[res$padj < 0.01 ]



RUVg.rs <- RUVg(x = round(2^expr.matrix-1),cIdx = ifelse(rownames(expr.matrix) %in% de.gene,FALSE,TRUE),k=1)
W <- RUVg.rs$W
rownames(W) <- colnames(expr.matrix)
df$W <- W[,1]
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition + W)

dds <- DESeq(dds)
res <- results(dds,contrast=c('condition','MET500','cell.line'))
res <- res[complete.cases(res),]
res <- res[order(res$pvalue),]


Basal.MET500.vs.CCLE.up.de.gene   <- (res)[res$padj < 0.001 & res$log2FoldChange >  1,]
Basal.MET500.vs.CCLE.down.de.gene <- (res)[res$padj < 0.001 & res$log2FoldChange < -1,]




############ DE analysis for non-Basal-organoid vs non-Basal-MET500 samples#########

expr.matrix  <- cbind(organoid.log2.read.count.matrix[common.gene,non.Basal.organoid],MET500.log2.read.count.matrix[common.gene,MET500.breast.cancer.polyA.non.Basal.sample])
df           <- data.frame(condition=c(rep(x='cell.line',times=length(non.Basal.organoid)), rep(x='MET500',times=length(MET500.breast.cancer.polyA.non.Basal.sample))))
df$condition <- factor(df$condition,levels = c('cell.line','MET500'))
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res <- res[complete.cases(res),]
res <- res[order(res$pvalue),]

de.gene   <- rownames(res)[res$padj < 0.01 ]



RUVg.rs <- RUVg(x = round(2^expr.matrix-1),cIdx = ifelse(rownames(expr.matrix) %in% de.gene,FALSE,TRUE),k=1)
W <- RUVg.rs$W
rownames(W) <- colnames(expr.matrix)
df$W <- W[,1]
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition + W)

dds <- DESeq(dds)
res <- results(dds,contrast=c('condition','MET500','cell.line'))
res <- res[complete.cases(res),]
res <- res[order(res$pvalue),]
non.Basal.MET500.vs.organoid.up.de.gene   <- (res)[res$padj < 0.001 & res$log2FoldChange >  1,]
non.Basal.MET500.vs.organoid.down.de.gene <- (res)[res$padj < 0.001 & res$log2FoldChange < -1,]


############ DE analysis for Basal-organoid vs Basal-MET500 samples#########

expr.matrix  <- cbind(organoid.log2.read.count.matrix[common.gene,Basal.organoid],MET500.log2.read.count.matrix[common.gene,MET500.breast.cancer.polyA.Basal.sample])
df           <- data.frame(condition=c(rep(x='cell.line',times=length(Basal.organoid)), rep(x='MET500',times=length(MET500.breast.cancer.polyA.Basal.sample))))
df$condition <- factor(df$condition,levels = c('cell.line','MET500'))
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res <- res[complete.cases(res),]
res <- res[order(res$pvalue),]

de.gene   <- rownames(res)[res$padj < 0.01 ]



RUVg.rs <- RUVg(x = round(2^expr.matrix-1),cIdx = ifelse(rownames(expr.matrix) %in% de.gene,FALSE,TRUE),k=1)
W <- RUVg.rs$W
rownames(W) <- colnames(expr.matrix)
df$W <- W[,1]
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition + W)

dds <- DESeq(dds)
res <- results(dds,contrast=c('condition','MET500','cell.line'))
res <- res[complete.cases(res),]
res <- res[order(res$pvalue),]
Basal.MET500.vs.organoid.up.de.gene   <- (res)[res$padj < 0.001 & res$log2FoldChange >  1,]
Basal.MET500.vs.organoid.down.de.gene <- (res)[res$padj < 0.001 & res$log2FoldChange < -1,]





########################
non.Basal.common.up.gene   <- intersect(non.Basal.MET500.vs.CCLE.up.de.gene %>% rownames,non.Basal.MET500.vs.organoid.up.de.gene %>% rownames) 
non.Basal.common.down.gene <- intersect(non.Basal.MET500.vs.CCLE.down.de.gene %>% rownames,non.Basal.MET500.vs.organoid.down.de.gene %>% rownames) 
Basal.common.up.gene       <- intersect(Basal.MET500.vs.CCLE.up.de.gene %>% rownames,  Basal.MET500.vs.organoid.up.de.gene %>% rownames) 
Basal.common.down.gene     <- intersect(Basal.MET500.vs.CCLE.down.de.gene %>% rownames,Basal.MET500.vs.organoid.down.de.gene %>% rownames) 
common.up.gene             <- intersect(non.Basal.common.up.gene,  Basal.common.up.gene)
common.down.gene           <- intersect(non.Basal.common.down.gene,Basal.common.down.gene)
