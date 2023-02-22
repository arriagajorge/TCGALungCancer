setwd("/home/jvasquez/Documents/TCGA-lung")
#FROM:
##  HANDS-ON:   NGS ANALYSIS  -->  Quantification files
## By Carlos Martínez Mira, Nov-2016
## By Sonia Tarazona, Nov-2016
##https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
## &&
## Data preparation: 
##      -Quality Control & bias removal
## By Cristobal Fresno - cristobalfresno@gmail.com
library(SummarizedExperiment)
library(TCGAbiolinks)
library(biomaRt)  
subtypeLUSC=read.table("subtypeLUSC.tsv",header=T, sep="\t")

xprssnLUSC <- GDCquery(project = "TCGA-LUSC",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "STAR - Counts",
                       barcode=subtypeLUSC$samples)

df <- xprssnLUSC[[1]][[1]]
GDCdownload(xprssnLUSC)
expreLUSCT <- GDCprepare(xprssnLUSC, summarizedExperiment = T)
expreLUSC <- GDCprepare(xprssnLUSC, summarizedExperiment = F)

tempLUSC = as.matrix(expreLUSC[,2:ncol(expreLUSC)])
rownames(tempLUSC) = expreLUSC$gene_id
expreLUSC = tempLUSC

stranded_firstLUSC <- expreLUSCT@assays@data@listData[["stranded_first"]]
unstrandedLUSC <- expreLUSCT@assays@data@listData[["unstranded"]]
stranded_secondLUSC <- expreLUSCT@assays@data@listData[["stranded_second"]]
tpm_unstrandLUSC <- tpm_unstrandLUSC <- expreLUSCT@assays@data@listData[["tpm_unstrand"]]
fpkm_unstrandLUSC <- expreLUSCT@assays@data@listData[["fpkm_unstrand"]]
fpkm_uq_unstrandLUSC <- expreLUSCT@assays@data@listData[["fpkm_uq_unstrand"]]
# "unstranded" "stranded_f" "stranded_s" "tpm_unstra" "fpkm_unstr" "fpkm_uq_un"

colnames(stranded_firstLUSC) <- expreLUSCT@colData@rownames
colnames(unstrandedLUSC) <- expreLUSCT@colData@rownames
colnames(stranded_secondLUSC) <- expreLUSCT@colData@rownames
colnames(tpm_unstrandLUSC) <- expreLUSCT@colData@rownames
colnames(fpkm_unstrandLUSC) <- expreLUSCT@colData@rownames
colnames(fpkm_uq_unstrandLUSC) <- expreLUSCT@colData@rownames

# subtype to duplicates #one only
i = substr(colnames(fpkm_unstrandLUSC), 1, 19)
j = i[duplicated(i)]
designExpLUSC=subtypeLUSC[c(which(!subtypeLUSC$samples%in%j),
                            as.numeric(sapply(which(subtypeLUSC$samples%in%j),rep,2))),]
dim(designExpLUSC); dim(fpkm_unstrandLUSC)
designExpLUSC=designExpLUSC[order(match(designExpLUSC$samples,substr(colnames(expreLUSC),1,19))),]
dim(designExpLUSC); dim(fpkm_unstrandLUSC)
designExpLUSC$barcode=colnames(fpkm_unstrandLUSC)

variables_ <- list(stranded_firstLUSC, unstrandedLUSC, stranded_secondLUSC, tpm_unstrandLUSC, 
                   fpkm_unstrandLUSC, fpkm_uq_unstrandLUSC)

write.table(fpkm_unstrandLUSC, "RNAseqLUSC.tsv", sep = '\t', quote = F)
#stranded_firstLUSC <- cbind(gene_id = expreLUSC$gene_id[1:60660], stranded_firstLUSC)
#stranded_firstLUSC <- cbind(gene_name = expreLUSC$gene_name[1:60660], stranded_firstLUSC)
#stranded_firstLUSC <- cbind(gene_type = expreLUSC$gene_type[1:60660], stranded_firstLUSC)


dim(designExpLUSC); dim(stranded_firstLUSC)

# expreLUSC[,"gene_type"][1:5]
# colnames(expreLUSC)[1:5]

add_gene_info <- function(dfLUSC){
  dfLUSC <- cbind(gene_type = expreLUSC[,"gene_type"][1:60660], dfLUSC)
  dfLUSC <- cbind(gene_name = expreLUSC[,"gene_name"][1:60660], dfLUSC)
  dfLUSC <- cbind(gene_id = rownames(expreLUSC)[1:60660], dfLUSC)
  return(dfLUSC)
}
unstrandedLUSC <- add_gene_info(unstrandedLUSC)
stranded_firstLUSC <- add_gene_info(stranded_firstLUSC)
stranded_secondLUSC <- add_gene_info(stranded_secondLUSC)
tpm_unstrandLUSC <- add_gene_info(tpm_unstrandLUSC)
fpkm_unstrandLUSC <- add_gene_info(fpkm_unstrandLUSC)
fpkm_uq_unstrandLUSC <- add_gene_info(fpkm_uq_unstrandLUSC)

dim(designExpLUSC)
# head(i); head(j); length(i); length(j); unique(j)

# keep only tenscript id not version numbers
rownames(unstrandedLUSC) <- unstrandedLUSC[,"gene_id"]
rownames(unstrandedLUSC) <- sapply(strsplit(rownames(unstrandedLUSC), ".", fixed=T),
                                   function(x) x[1])

rownames(stranded_firstLUSC) <- rownames(unstrandedLUSC)
rownames(stranded_secondLUSC) <- rownames(unstrandedLUSC)
rownames(tpm_unstrandLUSC) <- rownames(unstrandedLUSC)
rownames(fpkm_unstrandLUSC) <- rownames(unstrandedLUSC)
rownames(fpkm_uq_unstrandLUSC) <- rownames(unstrandedLUSC)

## "unstranded" "stranded_f" "stranded_s" "tpm_unstra" "fpkm_unstr" "fpkm_uq_un"
#annnotate GC content, length & biotype per transcript
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", 
                             "percentage_gene_gc_content", "gene_biotype",
                             "start_position","end_position","hgnc_id","hgnc_symbol"),
              filters = "ensembl_gene_id", 
              values=rownames(fpkm_unstrandLUSC),mart=mart) #its valid for every variable
#que son los espacios en blanco en my annot??
myannot$length=abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation
myannot=myannot[myannot$gene_biotype=="protein_coding"&
                  myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
exprots_hgnc=fpkm_unstrandLUSC[rownames(fpkm_unstrandLUSC)%in%myannot$ensembl_gene_id,]
exprots_hgnc <- exprots_hgnc[,4:ncol(exprots_hgnc)]
dim(exprots_hgnc)
#exprots_hgnc[,"gene_id"]
#[1] 19400   75

##check duplicated probes                           xd????
#myannot[myannot$hgnc_id == myannot$hgnc_id[duplicated(myannot$hgnc_id)],]
myannot2 <- myannot[unique(rownames(myannot)),]
dim(myannot2); dim(myannot)
# > myannot$hgnc_id[duplicated(myannot$hgnc_id)]
# [1] "HGNC:30046" "HGNC:11582" "HGNC:33853" "HGNC:4876" 
which(myannot2$hgnc_id == "HGNC:30046"); which(myannot2$hgnc_id == "HGNC:11582")
which(myannot2$hgnc_id == "HGNC:33853"); which(myannot2$hgnc_id == "HGNC:4876")

myannot2[c(18566,18728),]; myannot2[c(19327,19336),]
myannot2[c(12080,19340),]; myannot2[c(7581,19375),]

myannot3 <- myannot2[-c(18566,19327,19340,7581),]
dim(myannot2); dim(myannot3)
# length(unique(rownames(stranded_firstLUSC)))
# length(unique(rownames(stranded_firstLUSC)))
# [1] 60616

fpkm_LUSC <- fpkm_unstrandLUSC[unique(rownames(stranded_firstLUSC)),]

##################CHECK BIASES########################################################
library(NOISeq)
library(edgeR)

expre_hgnc

#format data for noiseq
noiseqData = readData(data = exprots_hgnc[unique(rownames(exprots_hgnc)),],
                      gc = myannot3[,1:2],
                      biotype = myannot3[,c(1,3)],factors=designExpLUSC,
                      length=myannot3[,c(1,8)])

noiseqData = readData(data = exprots_hgnc[unique(rownames(exprots_hgnc)),],
                      gc = myannot[,1:2],
                      biotype = myannot[,c(1,3)],factors=designExpLUSC,
                      length=myannot[,c(1,8)])

#1)check expression bias per subtype
mycountsbio = dat(noiseqData, type = "countsbio", factor = "subtype")
mycountsbio2 = dat(noiseqData, type = "countsbio")
mycountsbio3 = dat(noiseqData, factor = "subtype")
mycountsbio4 = dat(noiseqData)
unique(designExpLUSC$subtype)

png("CountsOri.png")