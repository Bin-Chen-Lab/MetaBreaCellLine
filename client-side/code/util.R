##### I put some functions  that will be used by other files##########

#####Function to pick out cell line##########
pick.out.cell.line <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
    marker.gene           <- intersect(rownames(expr.of.samples),(marker.gene))  
    marker.gene           <- intersect(rownames(expr.of.cell.lines),(marker.gene)) 
    correlation.matrix    <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
    cell.line.median.cor  <- apply(correlation.matrix,2,median) %>% sort(decreasing = TRUE)
    best.cell.line        <- names(cell.line.median.cor)[1]
    p.value.vec           <- foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% {
        v                     <- correlation.matrix[,cell.line]
        p.value               <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
    }
    names(p.value.vec) <- setdiff(names(cell.line.median.cor),best.cell.line)
    fdr.vec            <- p.adjust(p.value.vec,method='fdr')
    list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix )
}



# ####### Compute CCLE rna seq marker genes ##########
# load('server-side/RData/CCLE.RData')
# CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
# CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
# tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
# tmp.rank                    <- apply(tmp,2,rank)
# rank.mean                   <- apply(tmp.rank,1,mean)
# rank.sd                     <- apply(tmp.rank,1,sd)
# plot(x=rank.mean,y=rank.sd)
# lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
# CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
# CCLE.rna.seq.marker.gene.2000                 <- names(sort(rank.sd,decreasing =TRUE))[1:2000]
# CCLE.rna.seq.marker.gene.3000                 <- names(sort(rank.sd,decreasing =TRUE))[1:3000]
# 
# 
# 
# ##########
