# @title getData_lungCancer TCGA-LUAD
# @date 02/2023
# @author jvasquez

setwd("/home/jvasquez/Documents/TCGA-lung")
library(SummarizedExperiment)
library(TCGAbiolinks)
library(VennDiagram)
library(data.table)
library(tidyverse)

mthyltnLUAD <-  GDCquery(project = "TCGA-LUAD",
                         data.category = "DNA Methylation",
                         platform="Illumina Human Methylation 450")

df1 <- mthyltnLUAD[[1]][[1]]
mthyltnLUAD = getResults(mthyltnLUAD)

i=substr(mthyltnLUAD$cases,1,19)

xprssnLUAD <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type="STAR - Counts")#
xprssnLUAD=getResults(xprssnLUAD)

j=substr(xprssnLUAD$cases,1,19)
mirnasLUAD <- GDCquery(project = "TCGA-LUAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "miRNA Expression Quantification")
mirnasLUAD=getResults(mirnasLUAD)
k=substr(mirnasLUAD$cases,1,19)

##############CONCOURRENT MEASURES########################
sapply(list(i,j,k),function(x) length(unique(x)))
# [1] 503 589 564
#only samples with concurrent measurements are useful
samples=intersect(intersect(i,j),k)
length(samples)
#[1] 467

samples_df=data.frame(cbind(samples,
                            sapply(samples,function(x) 
                              unique(as.character(mthyltnLUAD$sample_type[i==x])))))
colnames(samples_df)[2]="tissue"
# unique(samples_df$tissue)
# [1] "Primary Tumor"       "Solid Tissue Normal" "Recurrent Tumor" 
samples_df$patient=substr(samples_df$samples,1,12) 

samples <- samples_df

subtypes=TCGAquery_subtype(tumor="LUAD")#subtype per patient
sum(samples$patient%in%subtypes$patient)
# [1] 193
# unique(subtypes$expression_subtype)
# [1] prox.-inflam prox.-prolif. TRU

#only classified samples  are useful
samples=samples[samples$patient%in%subtypes$patient,]
table(subtypes$expression_subtype[subtypes$patient%in%samples$patient])
# prox.-inflam prox.-prolif.           TRU 
#           68            48            70 

#normal subtype is discarded NOT normal tissue 
samples_2<-samples[which(!samples$patient%in%subtypes$patient[   #samples_2 is only a name for dont overwrite
  subtypes$expression_subtype=="Normal"]|
    samples$tissue=="Solid Tissue Normal"),] #there not are Subtype Normal

temp=table(samples_2[samples_2$patient%in%samples_2$patient[duplicated(
  samples_2$patient)],2:3])
table(apply(temp,2,paste,collapse=""))

#subtype per sample
samples_2$subtype=sapply(samples_2$patient,function(x) 
  subtypes$expression_subtype[subtypes$patient==x])

levels(samples_2$subtype) <- c(levels(samples_2$subtype), "normal")
samples_2 <- samples_2 %>% 
  mutate(subtype = replace(subtype, tissue =="Solid Tissue Normal", "normal"))

#samples_2["subtype"][samples_2["tissue"]=="Solid Tissue Normal"]="Normal"

table(samples_2$subtype)
#prox.-inflam prox.-prolif.           TRU        normal 
#         69            48            71             5 

write.table(samples_2,"subtypeLUAD.tsv",sep='\t',quote=F,row.names=F)

getwd()
#plot intersections
i=i[mthyltnLUAD$sample_type!="Solid Tissue Normal"]
j=j[xprssnLUAD$sample_type!="Solid Tissue Normal"]
k=k[mirnasLUAD$sample_type!="Solid Tissue Normal"]
levels_lung <- as.character(unique(subtypes$expression_subtype))

lista <- lapply(levels_lung,function(x)
  venn.diagram(x=list(
    A=unique(j[substr(j,1,12)%in%subtypes$patient[subtypes$expression_subtype==x]]),
    B=unique(i[substr(i,1,12)%in%subtypes$patient[subtypes$expression_subtype==x]]),
    C=unique(k[substr(k,1,12)%in%subtypes$patient[subtypes$expression_subtype==x]])),
    col = "transparent", fill = c("cornflowerblue","green","red"),
    alpha = 0.50,cex = 1.5,cat.cex = 1.5,margin = 0.1,
    fontfamily = "sans",cat.fontfamily=rep("sans",3),
    category.names=c("RNAseq","HM450","miRNAseq"),filename=x)) #estos nombres estan bien?

venn.diagram(x=list(
  A=substr(mthyltnLUAD$cases[mthyltnLUAD$sample_type=="Solid Tissue Normal"],1,19),
  B=substr(xprssnLUAD$cases[xprssnLUAD$sample_type=="Solid Tissue Normal"],1,19),
  C=substr(mirnasLUAD$cases[mirnasLUAD$sample_type=="Solid Tissue Normal"],1,19)),
  col = "transparent", fill = c("cornflowerblue","green","red"),
  alpha = 0.50,cex = 1.5,cat.cex = 1.5,margin = 0.1,
  fontfamily = "sans",cat.fontfamily=rep("sans",3),
  category.names=c("RNAseq","HM450","miRNAseq"),filename="Normal")


###############ADD CLINICAL INFO##########################
clin <- GDCquery_clinic("TCGA-LUAD","clinical")
#change tumor_stage by ajcc_pathologic_stage #jvasquez
clin <- clin[,c("bcr_patient_barcode","gender",
                "ajcc_pathologic_stage","race","vital_status")]
samples_2=cbind(samples_2,t(sapply(samples_2$patient,function(x) 
  clin[clin$bcr_patient_barcode==x,2:4])))
table(clin$gender) # change subtype by clin
#female   male 
#   280    242
table(clin$ajcc_pathologic)
#    Stage I   Stage IA   Stage IB   Stage II  Stage IIA  Stage IIB Stage IIIA Stage IIIB 
#         5        134        140          1         50         73         74         11 
# Stage IV 
# 26
table(clin$race)
# american indian or alaska native                            asian 
#                                 1                                8 
# black or african american                     not reported 
#                       53                               67 
# white 
# 393 
samples <- samples_2