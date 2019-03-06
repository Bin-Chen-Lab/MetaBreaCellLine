### R code to determine MET500 breast cancer subtypes ###


require(dplyr)
require(genefu)
load('server-side/RData/MET500.RData')


flag                                  <- grepl(x=MET500.sample.meta$cancer.type,pattern='Breast Invasive Ductal Carcinoma') & MET500.sample.meta$LibrarySelection == 'PolyA'
MET500.breast.cancer.polyA.sample     <- rownames(MET500.sample.meta)[flag] %>% as.character
flag                                  <- grepl(x=MET500.sample.meta$cancer.type,pattern='Breast Invasive Ductal Carcinoma') & MET500.sample.meta$LibrarySelection == 'Hybrid Selection'
MET500.breast.cancer.hybrid.sample    <- rownames(MET500.sample.meta)[flag] %>% as.character


#### Perform pam50 subtyping for MET500 polyA samples###########
pam50.gene.df               <- read.delim("~/Project/Cancer2CellLine/client-side/meta.data/pam50_entrez_gene_id_and_ensemble_gene_id.txt", stringsAsFactors=FALSE,comment.char = '#')
colnames(pam50.gene.df)[1]  <- 'EntrezGene.ID'
colnames(pam50.gene.df)[2]  <- 'probe' #damn it, genefu package doesnot mentioning this!
pam50.gene.df$EntrezGene.ID <- as.character(pam50.gene.df$EntrezGene.ID)
pam50.gene.expr             <- MET500.log2.fpkm.matrix[pam50.gene.df$probe %>% as.character,MET500.breast.cancer.polyA.sample] %>% t 
annot.matrix                <- pam50.gene.df[,1:2] %>% as.matrix
rownames(annot.matrix)      <- annot.matrix[,'probe']
pam50.subtype.rs            <- intrinsic.cluster.predict(sbt.model = pam50.robust,data = pam50.gene.expr,annot = annot.matrix,mapping=annot.matrix ,do.mapping = TRUE )

MET500.breast.cancer.polyA.LumB.sample    <- MET500.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'LumB']
MET500.breast.cancer.polyA.Basal.sample   <- MET500.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Basal']
MET500.breast.cancer.polyA.Her2.sample    <- MET500.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Her2']
MET500.breast.cancer.polyA.LumA.sample    <- MET500.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'LumA']
MET500.breast.cancer.polyA.Normal.sample  <- MET500.breast.cancer.polyA.sample[pam50.subtype.rs$subtype == 'Normal']

save(file = 'client-side/output/MET500.breast.cancer.meta.R.output/MET500.breast.cancer.meta.RData',
     list = c('MET500.breast.cancer.polyA.LumB.sample',  'MET500.breast.cancer.polyA.LumA.sample',
              'MET500.breast.cancer.polyA.Her2.sample',  'MET500.breast.cancer.polyA.Basal.sample',
              'MET500.breast.cancer.polyA.Normal.sample','MET500.breast.cancer.polyA.sample','MET500.breast.cancer.hybrid.sample'
             )
)

