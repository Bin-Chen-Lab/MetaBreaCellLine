
load('server-side/RData/CCLE.RData')
load('server-side/RData/microarray.RData')

# The gene copy number file is ONLY used to determine which CCLE cell lines were genotyped, I will NOT use the data in the file
CCLE_copynumber_byGene_2013.12.03 <- read.delim("~/Project/Cancer2CellLine/client-side/Data/CNV/CCLE_copynumber_byGene_2013-12-03.txt", stringsAsFactors=FALSE)
CCLE.genotyped.cell.line          <- colnames(CCLE_copynumber_byGene_2013.12.03)[6:ncol(CCLE_copynumber_byGene_2013.12.03)]
CCLE.genotyped.cell.line          <- gsub(x=CCLE.genotyped.cell.line,pattern = "^X",replacement = '')

# CCLE cell lines whole expression was profiled by rna-seq
CCLE.rna.seq.cell.line            <- colnames(CCLE.log2.rpkm.matrix)

# We only care about cell lines which have both genotype and expression data 
tmp                               <- intersect(CCLE.genotyped.cell.line,CCLE.rna.seq.cell.line)
CCLE.breast.cancer.cell.line      <- tmp[grepl(x=tmp,pattern='BREAST')]
CCLE.genotyped.and.rna.seqqed.cell.line                    <- tmp

# Divide breast cancer cell lines into metastatic and non-metastatic groups
CCLE.breast.cancer.cell.line.characteristic                <- read.csv("~/Project/Cancer2CellLine/client-side/meta.data/CCLE.breast.cancer.cell.line.characteristic.csv", stringsAsFactors=FALSE)
CCLE.breast.cancer.cell.line.characteristic$Cell.line.name <- paste0(CCLE.breast.cancer.cell.line.characteristic$Cell.line.name,sep='_BREAST')
flag                                                       <- CCLE.breast.cancer.cell.line.characteristic$Cell.line.name %in% colnames(CCLE.log2.rpkm.matrix)
CCLE.breast.cancer.cell.line.characteristic                <- CCLE.breast.cancer.cell.line.characteristic[flag,]
rownames(CCLE.breast.cancer.cell.line.characteristic)      <- CCLE.breast.cancer.cell.line.characteristic$Cell.line.name
CCLE.breast.cancer.metastatic.cell.line                    <- CCLE.breast.cancer.cell.line.characteristic$Cell.line.name[CCLE.breast.cancer.cell.line.characteristic$is.metastatic == 'yes']
CCLE.breast.cancer.non.metastatic.cell.line                <- CCLE.breast.cancer.cell.line.characteristic$Cell.line.name[CCLE.breast.cancer.cell.line.characteristic$is.metastatic == 'no']



save(file='client-side/output/CCLE.breast.cancer.cell.line.meta.R.output/CCLE.breast.cancer.cell.line.meta.RData',
     list=c('CCLE.breast.cancer.cell.line','CCLE.breast.cancer.metastatic.cell.line','CCLE.breast.cancer.non.metastatic.cell.line','CCLE.genotyped.and.rna.seqqed.cell.line')
     )