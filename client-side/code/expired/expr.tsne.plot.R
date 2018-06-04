load('client-side/output/combine.expression.R.output/combined.expr.matrix.RData')
load('client-side/output/combine.expression.R.output/normalized.combined.expr.matrix.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')
load('client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData')
require(foreach)
require(dplyr)



data.source                                                           <- rep('CCLE',times=ncol(combined.expr.matrix))
data.source[grepl(x=colnames(combined.expr.matrix),pattern='SRR')]    <- 'MET500'
data.source[grepl(x=colnames(combined.expr.matrix),pattern='TCGA')]   <- 'TCGA'
data.source[grepl(x=colnames(combined.expr.matrix),pattern='BREAST')] <- 'CCLE.breast.cancer'



normalized.data <- combat.normalized.combined.expr.matrix
# normalized.data <- peer.5.normalized.combined.expr.matrix
# normalized.data <- peer.3.normalized.combined.expr.matrix
normalized.data <- peer.normalized.combined.expr.matrix
normalized.data <- pca.normalized.combined.expr.matrix
#normalized.data <- combined.expr.matrix


CCLE.log2.rpkm.matrix                         <- normalized.data[,CCLE.cell.line] 
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
CCLE.rna.seq.marker.gene.2000                 <- names(sort(rank.sd,decreasing =TRUE))[1:2000]
CCLE.rna.seq.marker.gene.3000                 <- names(sort(rank.sd,decreasing =TRUE))[1:3000]
dist.obj <- 1- cor(normalized.data[CCLE.rna.seq.marker.gene.1000,] ,method='spearman')



CCLE.log2.rpkm.matrix                         <- normalized.data[,CCLE.cell.line] 
CCLE.var                                      <- apply(CCLE.log2.rpkm.matrix,1,mad) 
CCLE.rna.seq.marker.gene.100                 <- names(sort(CCLE.var,decreasing =TRUE))[1:100]

CCLE.rna.seq.marker.gene.1000                 <- names(sort(CCLE.var,decreasing =TRUE))[1:1000]
CCLE.rna.seq.marker.gene.2000                 <- names(sort(CCLE.var,decreasing =TRUE))[1:2000]
CCLE.rna.seq.marker.gene.3000                 <- names(sort(CCLE.var,decreasing =TRUE))[1:3000]
CCLE.rna.seq.marker.gene.5000                 <- names(sort(CCLE.var,decreasing =TRUE))[1:5000]
dist.obj <- dist(normalized.data[CCLE.rna.seq.marker.gene.100,] %>% t)




source('client-side/code/util.R')
require(foreach)
rs <- pick.out.cell.line(expr.of.cell.lines = normalized.data [,CCLE.cell.line],expr.of.samples = normalized.data [,data.source == 'MET500'],marker.gene = CCLE.rna.seq.marker.gene.1000)



require(Rtsne)
tsne.rs <- Rtsne(dist.obj)
# tsne.rs <- Rtsne(normalized.data[CCLE.rna.seq.marker.gene.1000,] %>% t)
draw.df           <- data.frame(dim1=tsne.rs$Y[,1],dim2=tsne.rs$Y[,2],data.source = data.source)
rownames(draw.df) <- colnames(normalized.data)

draw.df$subtype <- 'none'
flag <- rownames(draw.df) %in% MET500.breast.cancer.polyA.Basal.sample
draw.df[flag,'subtype'] <- 'Basal'
flag <- rownames(draw.df) %in% MET500.breast.cancer.polyA.Her2.sample
draw.df[flag,'subtype'] <- 'Her2'
flag <- rownames(draw.df) %in% MET500.breast.cancer.polyA.LumA.sample
draw.df[flag,'subtype'] <- 'LumA'
flag <- rownames(draw.df) %in% MET500.breast.cancer.polyA.LumB.sample
draw.df[flag,'subtype'] <- 'LumB'

draw.df[names(TCGA.breast.cancer.pam50.subtype),'subtype'] <- TCGA.breast.cancer.pam50.subtype
ggplot(draw.df) + geom_point(aes(x=dim1,y=dim2,shape=data.source,color=subtype),size=3) + theme_linedraw(base_size = 35,base_family = 'Arial')



# pca.rs <- prcomp(combined.expr.matrix %>% t)
# 
# 
# 
# 
# s <- (revmap(org.Hs.egENSEMBL) %>% as.list)
# gg <- s[CCLE.rna.seq.marker.gene.1000] %>% unlist
# 
# 
# common.genes <- intersect(CCLE.rna.seq.marker.gene.1000,rownames(combine.expr.matrix))
# dist.obj <- (1- cor(combine.expr.matrix[common.genes,c(MET500.breast.cancer.polyA.sample,CCLE.breast.cancer.cell.line,colnames(TCGA.breast.cancer.log2.fpkm.matrix) )],method='spearman')) %>% as.dist
# common.genes <- intersect(gg,rownames(normalized.combine.expr.matrix))
# 
# dist.obj <- (1- cor(normalized.combine.expr.matrix[common.genes,c(MET500.breast.cancer.polyA.sample,CCLE.breast.cancer.cell.line )],method='spearman')) %>% as.dist
# 
# 
# 
# col.annotataion <- data.frame( data.source = c(rep(x='MET500',times=length(MET500.breast.cancer.polyA.sample)),
#                                                rep(x='CCLE',  times=ncol(CCLE.log2.rpkm.matrix) ),
#                                                rep(x='TCGA',  times=ncol(TCGA.breast.cancer.log2.fpkm.matrix))
# )
# )
# rownames(col.annotataion) <- colnames(combine.expr.matrix )
# 
# 
# 
# 
# tsne.rs  <- Rtsne(dist.obj)
# 
# rownames(tsne.rs$Y) <- colnames(normalized.combine.expr.matrix)
# rownames(tsne.rs$Y) <- c(MET500.breast.cancer.polyA.sample,CCLE.breast.cancer.cell.line,colnames(TCGA.breast.cancer.log2.fpkm.matrix))
# rownames(tsne.rs$Y) <- c(MET500.breast.cancer.polyA.sample,CCLE.breast.cancer.cell.line)
# 
# draw.df             <- data.frame(x=tsne.rs$Y[,1],y=tsne.rs$Y[,2],col=col.annotataion[rownames(tsne.rs$Y),'data.source'])
# rownames(draw.df)   <- rownames(tsne.rs$Y)
# draw.df$col         <- as.character(draw.df$col)
# draw.df[draw.df$col == 'CCLE','col']                      <- 'CCLE.non.breast.cancer'
# draw.df[grepl(x=rownames(draw.df),pattern='BREAST'),'col'] <- 'CCLE.breast.cancer'
# ggplot(draw.df) + geom_point(aes(x=x,y=y,col=col),size=3) + theme_gray(base_size = 25) + xlab('Dim1') + ylab('Dim2')

