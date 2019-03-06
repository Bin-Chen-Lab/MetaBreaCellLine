require(plyr)
require(dplyr)
require(foreach)
require(ggplot2)



# Organize organoid RNASeq data
BreastCancerBioBank_RNAseq_readCounts_RPKM           <- read.delim("~/Project/Cancer2CellLine/client-side/Data/organoid.expr/BCBB_RNAseq_SummaryLevelData/BreastCancerBioBank_RNAseq_readCounts_RPKM.txt", stringsAsFactors=FALSE)
rownames(BreastCancerBioBank_RNAseq_readCounts_RPKM) <- BreastCancerBioBank_RNAseq_readCounts_RPKM$X
BreastCancerBioBank_RNAseq_readCounts_RPKM$X         <- NULL
organoid.rpkm.matrix                                 <- as.matrix(BreastCancerBioBank_RNAseq_readCounts_RPKM)


BreastCancerBioBank_RNAseq_readCounts_raw             <- read.delim("~/Project/Cancer2CellLine/client-side/Data/organoid.expr/BCBB_RNAseq_SummaryLevelData/BreastCancerBioBank_RNAseq_readCounts_raw.txt", stringsAsFactors=FALSE)
rownames(BreastCancerBioBank_RNAseq_readCounts_raw)   <- BreastCancerBioBank_RNAseq_readCounts_raw$gene
BreastCancerBioBank_RNAseq_readCounts_raw$gene        <- NULL
organoid.read.count.matrix                            <- as.matrix(BreastCancerBioBank_RNAseq_readCounts_raw)



remove.version <- function(x) {
    tmp <- strsplit(x=x,split = "\\.") %>% unlist 
    tmp[1]
}
organoid.line                         <- sapply(colnames(organoid.rpkm.matrix),remove.version)
organoid.line[25]                     <- 'W894.10B01A0' #hmm, two replicates for W894?
organoid.line                         <- paste(organoid.line,'_BREAST_ORGANOID',sep='')
colnames(organoid.rpkm.matrix)        <- organoid.line
colnames(organoid.read.count.matrix)  <- organoid.line

organoid.log2.rpkm.matrix             <- log2(organoid.rpkm.matrix+1)
organoid.log2.read.count.matrix       <- log2(organoid.read.count.matrix+1)


require(genefu)
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- organoid.rpkm.matrix[pam50.gene.df$probe %>% as.character,] %>% t 
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )

LumA.organoid               <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'LumA']
LumB.organoid               <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'LumB']
Her2.organoid               <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'Her2']
Basal.organoid              <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'Basal']
Normal.organoid             <- colnames(organoid.rpkm.matrix)[pam50.subtype.rs$subtype == 'Normal']


save(file='client-side/output/organize.breast.cancer.organoid.data.R.output/organoid.RData',list=c('organoid.log2.rpkm.matrix','organoid.log2.read.count.matrix','LumA.organoid','LumB.organoid','Her2.organoid','Basal.organoid'))