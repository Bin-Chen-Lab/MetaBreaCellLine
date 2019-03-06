### This R code file performs gene expression correlation analysi for in-vivo MDAMB231.####
### Basicaly, we downloaded in vivo lung-metastasis-derived MDAMB231 dataset from GEO and merged them with CCLE breast cancer cell lines, then 
### we matched these cell lines with lung-derived metastatic breast caner samples

require(GEOquery)
require(foreach)
require(dplyr)
require(plyr)
require(reshape2)
load('client-side/output/correlation.analysis.with.expression.microarray.R.output/correlation.analysis.with.expression.microarray.RData')
load('client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData')



GSE2603.series.data           <- read.delim("client-side/Data/GSE2603_series_matrix.txt", stringsAsFactors=FALSE)
rownames(GSE2603.series.data) <- GSE2603.series.data$ID_REF
GSE2603.series.data$ID_REF    <- NULL
MDAMB231.microarray.expr.matrix      <- as.matrix(GSE2603.series.data)[,1:21]
cancer.sample.microarray.expr.matrix <- as.matrix(GSE2603.series.data)[,23:121]


load('server-side/RData/microarray.RData')
source('client-side/code/util.R')


### Compute 1000 most-varied probesets ###
tmp                         <- CCLE.microarray.expr.matrix
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.microarray.marker.gene.1000 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]


### Match cell lines to lung-derived metastatic breast cancer samples ###
combined.data     <-  cbind(CCLE.microarray.expr.matrix[,CCLE.breast.cancer.cell.line],MDAMB231.microarray.expr.matrix[CCLE.microarray.expr.matrix %>% rownames,])
biopsy.site       <-  metastatic.breast.cancer.microarray.meta.df$biopsy.site
lung.samples      <-  microarray.vs.CCLE.lung.correlation.result$correlation.matrix %>% rownames
basal.samples     <-  microarray.vs.CCLE.Basal.correlation.result$correlation.matrix %>% rownames


microarray.lung.vs.CCLE.plus.lung.metastasis.MDAMB231 <- pick.out.cell.line(expr.of.samples = metastatic.breast.cancer.microarray.expr.matrix[,(lung.samples)],
                                                                            expr.of.cell.lines = combined.data,
                                                                            marker.gene = CCLE.microarray.marker.gene.1000
                                                                           )

lung.derived.MDAMB231 <-  c('GSM49955','GSM49956','GSM49957','GSM50017','GSM50018','GSM50019','GSM50020','GSM50021','GSM50022')

save(file = 'client-side/output/double.check.MDAMB231.R.output/double.check.MDAMB231.RData',list=c('microarray.lung.vs.CCLE.plus.lung.metastasis.MDAMB231','lung.derived.MDAMB231'))






                          