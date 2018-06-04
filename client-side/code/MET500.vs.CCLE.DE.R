### R code to perform differential expression analysis between MET500 breast cancer samples and CCLE cell lines ###
### Note: in this DE analysis, the comparision was done between  LumA/LumB/Her2 MET500 samples and GOOD cell lines (see line 12)
### For basal-like MET500 samples, since we DID NOT detect good cell lines, the comparision was NOT done!


load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
load('client-side/output/pick.out.cell.line.by.subtype.R.output/pick.out.cell.line.by.subtype.RData')
load('client-side/output/MET500.breast.cancer.meta.R.output//MET500.breast.cancer.meta.RData')


MET500.breast.cancer.polyA.non.Basal.sample           <- c(MET500.breast.cancer.polyA.LumB.sample,MET500.breast.cancer.polyA.LumA.sample,MET500.breast.cancer.polyA.Her2.sample)
MET500.breast.cancer.polyA.Basal.sample               <- c(MET500.breast.cancer.polyA.Basal.sample)
non.Basal.good.cell.line                              <- c(MET500.LumA.good.cell.line,MET500.LumB.good.cell.line,MET500.Her2.good.cell.line) %>% unique

############ DE analysis for non-Basal-cell-lines vs non-basal-MET500 samples#########
require(DESeq2)
common.gene <- intersect(rownames(CCLE.log2.rpkm.matrix),rownames(MET500.log2.fpkm.matrix))
expr.matrix <- cbind(CCLE.log2.read.count.matrix[common.gene,non.Basal.good.cell.line],MET500.log2.read.count.matrix[common.gene,MET500.breast.cancer.polyA.non.Basal.sample])
flag        <- apply(expr.matrix,1,function(x) sum(x>=1)) # filter lowly expressed genes
expr.matrix <- expr.matrix[flag == ncol(expr.matrix),]


df <- data.frame(condition=c(rep(x='cell.line',times=length(non.Basal.good.cell.line)), rep(x='MET500',times=length(MET500.breast.cancer.polyA.non.Basal.sample))))
df$condition <- factor(df$condition,levels = c('cell.line','MET500'))
dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
                              colData = df,
                              design = ~ condition)

dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$pvalue),]
non.Basal.up.de.gene   <- (res)[res$padj < 0.001 & res$log2FoldChange >  1,]
non.Basal.down.de.gene <- (res)[res$padj < 0.001 & res$log2FoldChange < -1,]

# ####### Since there are NOT many availavle good Basal cell lines, we did not do DE analysis here##########
# expr.matrix <- cbind(CCLE.log2.read.count.matrix[common.gene,Basal.good.cell.line],MET500.log2.read.count.matrix[common.gene,MET500.breast.cancer.polyA.Basal.sample])
# flag        <- apply(expr.matrix,1,function(x) sum(x>=1)) # filter lowly expressed genes
# expr.matrix <- expr.matrix[flag == ncol(expr.matrix),]
# 
# 
# df <- data.frame(condition=c(rep(x='cell.line',times=length(Basal.good.cell.line)), rep(x='MET500',times=length(MET500.breast.cancer.polyA.Basal.sample))))
# df$condition <- factor(df$condition,levels = c('cell.line','MET500'))
# dds <- DESeqDataSetFromMatrix(countData = round(2^expr.matrix-1),
#                               colData = df,
#                               design = ~ condition)
# 
# dds <- DESeq(dds)
# res <- results(dds)
# res <- res[order(res$pvalue),]
# Basal.up.de.gene   <- (res)[res$padj < 0.001 & res$log2FoldChange >  1,]
# Basal.down.de.gene <- (res)[res$padj < 0.001 & res$log2FoldChange < -1,]


### Save the results ###
save(file='client-side/output/MET500.vs.CCLE.DE.R.output/MET500.vs.CCLE.DE.RData',list=c('non.Basal.up.de.gene','non.Basal.down.de.gene'))
write.csv(x=non.Basal.up.de.gene,  file='client-side/output/MET500.vs.CCLE.DE.R.output/non.Basal.up.de.gene.csv',  row.names=TRUE,quote=F)
write.csv(x=non.Basal.down.de.gene,file='client-side/output/MET500.vs.CCLE.DE.R.output/non.Basal.down.de.gene.csv',row.names=TRUE,quote=F)
