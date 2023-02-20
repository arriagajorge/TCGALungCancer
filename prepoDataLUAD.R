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
subtypeLUAD=read.table("subtypeLUAD.tsv",header=T, sep="\t")

xprssnLUAD <- GDCquery(project = "TCGA-LUAD",
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "STAR - Counts",
                   barcode=subtypeLUAD$samples)

df <- xprssnLUAD[[1]][[1]]
GDCdownload(xprssnLUAD)
expreLUADT <- GDCprepare(xprssnLUAD, summarizedExperiment = T)
expreLUAD <- GDCprepare(xprssnLUAD, summarizedExperiment = F)

tempLUAD = as.matrix(expreLUAD[,2:ncol(expreLUAD)])
rownames(tempLUAD) = expreLUAD$gene_id
expreLUAD = tempLUAD

write.table(expreLUAD, "RNAseqLUAD.tsv", sep = '\t', quote = F)

stranded_firstLUAD <- expreLUADT@assays@data@listData[["stranded_first"]]
unstrandedLUAD <- expreLUADT@assays@data@listData[["unstranded"]]
stranded_secondLUAD <- expreLUADT@assays@data@listData[["stranded_second"]]
tpm_unstrandLUAD <- tpm_unstrandLUAD <- expreLUADT@assays@data@listData[["tpm_unstrand"]]
fpkm_unstrandLUAD <- expreLUADT@assays@data@listData[["fpkm_unstrand"]]
fpkm_uq_unstrandLUAD <- expreLUADT@assays@data@listData[["fpkm_uq_unstrand"]]
# "unstranded" "stranded_f" "stranded_s" "tpm_unstra" "fpkm_unstr" "fpkm_uq_un"

variables_ <- list(stranded_firstLUAD, unstrandedLUAD, stranded_secondLUAD, tpm_unstrandLUAD, 
                fpkm_unstrandLUAD, fpkm_uq_unstrandLUAD)

colnames(stranded_firstLUAD) <- expreLUADT@colData@rownames
colnames(unstrandedLUAD) <- expreLUADT@colData@rownames
colnames(stranded_secondLUAD) <- expreLUADT@colData@rownames
colnames(tpm_unstrandLUAD) <- expreLUADT@colData@rownames
colnames(fpkm_unstrandLUAD) <- expreLUADT@colData@rownames
colnames(fpkm_uq_unstrandLUAD) <- expreLUADT@colData@rownames

#stranded_firstLUAD <- cbind(gene_id = expreLUAD$gene_id[1:60660], stranded_firstLUAD)
#stranded_firstLUAD <- cbind(gene_name = expreLUAD$gene_name[1:60660], stranded_firstLUAD)
#stranded_firstLUAD <- cbind(gene_type = expreLUAD$gene_type[1:60660], stranded_firstLUAD)

# subtype to duplicates #one only
i = substr(colnames(stranded_firstLUAD), 1, 19)
j = i[duplicated(i)]
designExpLUAD=subtypeLUAD[c(which(!subtypeLUAD$samples%in%j),
                    as.numeric(sapply(which(subtypeLUAD$samples%in%j),rep,2))),]
designExpLUAD=designExpLUAD[order(match(designExpLUAD$samples,substr(colnames(expreLUAD),1,19))),]
designExpLUAD$barcode=colnames(stranded_firstLUAD)

# expreLUAD[,"gene_type"][1:5]
# colnames(expreLUAD)[1:5]

add_gene_info <- function(dfLUAD){
  dfLUAD <- cbind(gene_type = expreLUAD[,"gene_type"][1:60660], dfLUAD)
  dfLUAD <- cbind(gene_name = expreLUAD[,"gene_name"][1:60660], dfLUAD)
  dfLUAD <- cbind(gene_id = rownames(expreLUAD)[1:60660], dfLUAD)
  return(dfLUAD)
}
unstrandedLUAD <- add_gene_info(unstrandedLUAD)
stranded_firstLUAD <- add_gene_info(stranded_firstLUAD)
stranded_secondLUAD <- add_gene_info(stranded_secondLUAD)
tpm_unstrandLUAD <- add_gene_info(tpm_unstrandLUAD)
fpkm_unstrandLUAD <- add_gene_info(fpkm_unstrandLUAD)
fpkm_uq_unstrandLUAD <- add_gene_info(fpkm_uq_unstrandLUAD)

dim(designExpLUAD)
# head(i); head(j); length(i); length(j); unique(j)

# keep only tenscript id not version numbers
rownames(unstrandedLUAD) <- unstrandedLUAD[,"gene_id"]
rownames(unstrandedLUAD) <- sapply(strsplit(rownames(unstrandedLUAD), ".", fixed=T),
                                   function(x) x[1])

rownames(stranded_firstLUAD) <- rownames(unstrandedLUAD)
rownames(stranded_secondLUAD) <- rownames(unstrandedLUAD)
rownames(tpm_unstrandLUAD) <- rownames(unstrandedLUAD)
rownames(fpkm_unstrandLUAD) <- rownames(unstrandedLUAD)
rownames(fpkm_uq_unstrandLUAD) <- rownames(unstrandedLUAD)

## "unstranded" "stranded_f" "stranded_s" "tpm_unstra" "fpkm_unstr" "fpkm_uq_un"
#annnotate GC content, length & biotype per transcript
mart=useEnsembl("ensembl",dataset="hsapiens_gene_ensembl")
myannot=getBM(attributes = c("ensembl_gene_id", 
                             "percentage_gene_gc_content", "gene_biotype",
                             "start_position","end_position","hgnc_id","hgnc_symbol"),
              filters = "ensembl_gene_id", 
              values=rownames(stranded_firstLUAD),mart=mart) #its valid for every variable
#que son los espacios en blanco en my annot??
myannot$length=abs(myannot$end_position-myannot$start_position)

#filter transcripts withouth annotation
myannot=myannot[myannot$gene_biotype=="protein_coding"&
                  myannot$hgnc_symbol!="",]
myannot=myannot[!duplicated(myannot$ensembl_gene_id),]
exprots_hgnc=stranded_firstLUAD[rownames(stranded_firstLUAD)%in%myannot$ensembl_gene_id,]
exprots_hgnc <- exprots_hgnc[,4:ncol(exprots_hgnc)]
dim(exprots_hgnc)
#exprots_hgnc[,"gene_id"]
#[1] 19400   201

##check duplicated probes                           xd????
#myannot[myannot$hgnc_id == myannot$hgnc_id[duplicated(myannot$hgnc_id)],]
myannot <- myannot[unique(rownames(myannot)),]
dim(myannot2); dim(myannot)
# > myannot$hgnc_id[duplicated(myannot$hgnc_id)]
# [1] "HGNC:30046" "HGNC:11582" "HGNC:33853" "HGNC:4876" 
which(myannot2$hgnc_id == "HGNC:30046"); which(myannot2$hgnc_id == "HGNC:11582")
which(myannot2$hgnc_id == "HGNC:33853"); which(myannot2$hgnc_id == "HGNC:4876")

myannot2[c(18566,18728),]; myannot2[c(19327,19336),]
myannot2[c(12080,19340),]; myannot2[c(7581,19375),]

myannot3 <- myannot2[-c(18566,19327,12080,7581),]
dim(myannot2); dim(myannot3)
# length(unique(rownames(stranded_firstLUAD)))
# length(unique(rownames(stranded_firstLUAD)))
# [1] 60616

str_first_LUAD <- stranded_firstLUAD[unique(rownames(stranded_firstLUAD)),]

##################CHECK BIASES########################################################
# Error in data.frame(numeric(n), row.names = nms) : 
#   duplicate row.names: ENSG00000002586, ENSG00000124333, ENSG00000124334, ENSG00000167393, ENSG00000169084, ENSG00000169093, ENSG00000169100, ENSG00000178605, ENSG00000182162, ENSG00000182378, ENSG00000182484, ENSG00000185291, ENSG00000185960, ENSG00000196433, ENSG00000197976, ENSG00000198223, ENSG00000205755, ENSG00000214717
# > length(unique(rownames(stranded_firstLUAD)))

library(NOISeq)
library(edgeR)

#format data for noiseq
noiseqData = readData(data = exprots_hgnc[unique(rownames(exprots_hgnc)),],
                      gc = myannot3[,1:2],
                      biotype = myannot3[,c(1,3)],factor=designExpLUAD,
                      length=myannot3[,c(1,8)])
#1)check expression bias per subtype
mycountsbio = dat(noiseqData, type = "countsbio", factor = "subtype")
#patients with repeated measures
png("CountsOri.png")
explo.plot(mycountsbio, plottype = "boxplot",samples = 1:5)
dev.off()
#2)check for low count genes
png("lowcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:5)
dev.off()
png("lowCountThres.png")
hist(rowMeans(cpm(exprots_hgnc,log=T)),ylab="genes",
     xlab="mean of log CPM",col="gray")
abline(v=0,col="red")
dev.off()