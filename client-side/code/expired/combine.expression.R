require(dplyr)
require(stringr)
load('server-side/RData/CCLE.RData')
load('server-side/RData/MET500.RData')
load('server-side/RData/TCGA.breast.cancer.RData')

flag                                  <- grepl(x=MET500.sample.meta$cancer.type,pattern='Breast Invasive Ductal Carcinoma') & MET500.sample.meta$LibrarySelection == 'PolyA'
MET500.breast.cancer.polyA.sample     <- rownames(MET500.sample.meta)[flag] %>% as.character

common.genes <- intersect(rownames(MET500.log2.fpkm.matrix),rownames(CCLE.log2.rpkm.matrix))
common.genes <- intersect(common.genes,rownames(TCGA.breast.cancer.log2.fpkm.matrix))

combined.expr.matrix <- cbind(MET500.log2.fpkm.matrix[common.genes,MET500.breast.cancer.polyA.sample],
                              CCLE.log2.rpkm.matrix[common.genes,],
                              TCGA.breast.cancer.log2.fpkm.matrix[common.genes,]
)



refseq.id.to.ensemble.gene           <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/refseq.id.to.ensemble.gene.txt", stringsAsFactors=FALSE)
colnames(refseq.id.to.ensemble.gene) <- c('ensemble.gene.id','ensemble.tr.id','refseq.id')
data                                 <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/human.house.keeping.gene.txt", header=FALSE, stringsAsFactors=FALSE)
refseq.id                            <- data$V2 %>% as.character %>% str_trim

flag                     <- match(refseq.id,refseq.id.to.ensemble.gene$refseq.id)
ensemble.gene.id         <- refseq.id.to.ensemble.gene$ensemble.gene.id[flag]
human.house.keeping.gene <- ensemble.gene.id[is.na(ensemble.gene.id) == FALSE]
human.house.keeping.gene <- human.house.keeping.gene[human.house.keeping.gene %in% common.genes]

save(file='client-side/output/combine.expression.R.output/combined.expr.matrix.RData',list=c('combined.expr.matrix','human.house.keeping.gene'))


