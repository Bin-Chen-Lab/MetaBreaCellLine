pam50.gene.df   <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')
pam50.gene     <- pam50.gene.df$ensemble.gene.id %>% as.character

cc <- intersect(CCLE.rna.seq.marker.gene.1000,rownames(MET500.log2.fpkm.matrix))

MET500.pam50.gene.expr             <- MET500.log2.fpkm.matrix[pam50.gene %>% as.character,MET500.breast.cancer.polyA.sample] 
CCLE.pam50.gene.expr               <- CCLE.log2.rpkm.matrix[pam50.gene %>% as.character,CCLE.breast.cancer.cell.line]

m <- cbind(MET500.pam50.gene.expr,CCLE.pam50.gene.expr)
dist.matrix <- 1-  cor(m,method='spearman')
pheatmap(dist.matrix,cluster_rows = T,cluster_cols = F)

CCLE.pam50.subtype <- pam50.subtype.rs.rna.seq$subtype
MET500.pam50.subtype <- pam50.subtype.rs$subtype

pam50.subtype <- c(CCLE.pam50.subtype,MET500.pam50.subtype)



require(tsne)
tsne.rs <- tsne(X = as.dist(dist.matrix))
#tsne.rs <- tsne(X = t(m))

rownames(tsne.rs) <- rownames(dist.matrix)

df <- data.frame(x=tsne.rs[,1],y=tsne.rs[,2],subtype=pam50.subtype[rownames(tsne.rs)],source=ifelse(grepl(x=rownames(tsne.rs),pattern='SRR'),'MET500','CCLE'))
ggplot(df,aes(x=x,y=y,col=subtype,shape=source)) + geom_point(size=5)+theme_gray(base_size = 30) + xlab('tsne-1') + ylab('tsne-2')



plot(tsne.rs,col=ifelse(rownames(tsne.rs) %in% MET500.breast.cancer.polyA.Basal.sample ,'red','black'),pch=19)
points(tsne.rs,col=ifelse(rownames(tsne.rs) %in% good.cell.line,'blue','black'),pch=19)
points(tsne.rs,col=ifelse(rownames(tsne.rs) %in% c('HCC70_BREAST','HCC2157_BREAST'),'purple','black'),pch=19)


plot(tsne.rs,col=ifelse(rownames(tsne.rs) %in% MET500.breast.cancer.polyA.sample,'red','black'),pch=19)
