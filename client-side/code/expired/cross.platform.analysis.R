require(plyr)
require(plyr)
load('client-side/output/correlation.analysis.with.expression.R.output/correlation.analysis.with.expression.RData')
load('client-side/output/correlation.analysis.with.expression.microarray.R.output/correlation.analysis.with.expression.microarray.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')


####### show analysis results with tissue specific data#######
plot(x=MET500.vs.CCLE.polyA.liver.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=MET500.vs.CCLE.polyA.lymph.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     xlab='Liver',ylab='Lymph node',cex=1.5,pch=19,cex.axis=1.5,cex.lab=1.5
    )
cor(x=MET500.vs.CCLE.polyA.liver.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=MET500.vs.CCLE.polyA.lymph.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     method='spearman'
)


############# show analysis results with subtype specific data #############
m <- cbind(MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
           MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
           MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
           MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
          )
           
colnames(m) <- c('LumA','LumB','Basal','Her2')
pairs(m,pch=19,cex=1.5,cex.axis=1.5,cex.lab=1.5)
cor.matrix <- cor(m,method='spearman')
pheatmap(cor.matrix,cluster_rows = F,cluster_cols = F)


########## correlation between subtype,microarray

m <-  cbind(microarray.vs.CCLE.LumA.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
            microarray.vs.CCLE.LumB.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
            microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
            microarray.vs.CCLE.Her2.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
)

colnames(m) <- c('LumA','LumB','Basal','Her2')
pairs(m,pch=19,cex=1.5,cex.axis=1.5,cex.lab=1.5)
cor.matrix <- cor(m,method='spearman')
pheatmap(cor.matrix,cluster_rows = F,cluster_cols = F)
############# correlation between tissue , microarray
m <- cbind(microarray.vs.CCLE.bone.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
            microarray.vs.CCLE.brain.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
            microarray.vs.CCLE.liver.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
            microarray.vs.CCLE.lung.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
            microarray.vs.CCLE.lymph.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line]
            
)

colnames(m) <- c('bone','brain','liver','lung','Lymph')
pairs(m,pch=19,cex=1.5,cex.axis=1.5,cex.lab=1.5)
cor.matrix <- cor(m,method='spearman')
pheatmap(cor.matrix,cluster_rows = F,cluster_cols = F)











###### compare two platforms, with ALL samples#######
plot(x=MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=microarray.vs.CCLE.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     xlab='MET500',ylab='microarray',xlim=c(0,0.5),ylim=c(0,0.5))
cor(x=MET500.vs.CCLE.polyA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
    y=microarray.vs.CCLE.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
    method='spearman')


###### compare two platforms, with liver samples#######
plot(x=MET500.vs.CCLE.polyA.liver.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=microarray.vs.CCLE.liver.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     xlab='MET500',ylab='microarray')
cor(x=MET500.vs.CCLE.polyA.liver.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=microarray.vs.CCLE.liver.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     method='spearman')



plot(x=MET500.vs.CCLE.polyA.lymph.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=microarray.vs.CCLE.lymph.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     xlab='MET500',ylab='microarray')
cor(x=MET500.vs.CCLE.polyA.lymph.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
    y=microarray.vs.CCLE.lymph.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
    method='spearman')



###### compare two platforms, with different subtypes #######
plot(x=MET500.vs.CCLE.polyA.LumA.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=microarray.vs.CCLE.LumA.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line])

plot(x=MET500.vs.CCLE.polyA.LumB.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=microarray.vs.CCLE.LumB.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line])

plot(x=MET500.vs.CCLE.polyA.Basal.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=microarray.vs.CCLE.Basal.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line])


plot(x=MET500.vs.CCLE.polyA.Her2.expression.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line],
     y=microarray.vs.CCLE.Her2.correlation.result$cell.line.median.cor[CCLE.breast.cancer.cell.line])





