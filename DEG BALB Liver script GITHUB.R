#----Packages-----
library(splitstackshape)
library(tidyr)
library(stringr)

library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(VennDiagram)
library(readxl)
library(calibrate)
library(HGNChelper)

library(ggplot2)
library(ggrepel)

#----Load RNAseq Data Set----
#Insert your directory here
count <- read.delim("C:/Users/jbeierle/Dropbox/lab/BALBc/BALBc QTL Mapping Fixed Genotyping/For GITHUB/Liver RNAseq/count_liver.txt")

count_working<-count_working[!grepl("Gm", count_working$id),]  #Remove unannotated genes
count_working<-count_working[!grepl("Rik", count_working$id),]

countdata<-count_working[,-c(2:6)] #Remove unneeded information

rownames(countdata) <- countdata[,1] #Make genes rownames
countdata <- countdata[,-(1)]

#-----QC to Remove lowly detected genes-----
myCPM <- cpm (countdata)
head(myCPM)

# >10 counts per mil for a library size of 30 mil on average
thresh <- myCPM > 0.30
head(thresh)
table(rowSums(thresh))
#keep genes that have at least 3 TRUES in each row of thresh - meaning genes that have >10 cpm in at least 3 samples
keep <- rowSums(thresh) >= 3
summary(keep)
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- countdata[keep,]
dim(counts.keep)

write.table(counts.keep, file="count_Liver_FinalQC.txt", sep = "\t", quote=FALSE,
            row.name=TRUE, col.names = TRUE)
#-----Subsetting By sex-----
counts.keep.F <- counts.keep[,c(1:8,17:19,21:25)]
sampleinfo.F <- sampleinfo[which(sampleinfo$Sex=="0"),]

counts.keep.M <- counts.keep[,c(9:16,20,26:32)]
sampleinfo.M <- sampleinfo[which(sampleinfo$Sex=="1"),]

.
#-------DEG analysis All mice------
sampleinfo <- read.delim("annotation_liver.txt")

y <- DGEList(counts.keep)
names(y)   # See what slots are stored in y
y$samples   # Library size information is stored in the samples slot


#----All Mice: EFFECT OF GENOTYPE, BLOCKING TX
y <- DGEList(counts = counts.keep)
y <- calcNormFactors(y)
Tx <- factor(sampleinfo$Tx)
Strain <- factor(sampleinfo$Strain)
design <- model.matrix(~Tx+Strain)

v <- voom(y,design,plot=TRUE)

fit <- lmFit(v, design)
fit_geno <- eBayes(fit)

Geno <- topTable(fit_geno,coef = "Strain1", n = "inf", sort = "p", adjust.method = "BH")

Geno$abs.logFC <- abs(Geno$logFC)
write.table(Geno, file="Geno-Tx_all J ref.txt", row.names = TRUE, sep="\t")
#At this point I have had to manually add the header "Gene" to the outputed text file
#file in col 1

#Insert your directory here
Geno.All <- read.delim("C:/Users/jbeierle/Dropbox/lab/BALBc/BALBc QTL Mapping Fixed Genotyping/For GITHUB/Liver RNAseq/Geno-Tx_all J ref.txt")


#-------DEG analysis FEMALE mice------
y <- DGEList(counts.keep.F)
names(y)   # See what slots are stored in y
y$samples   # Library size information is stored in the samples slot


#----EFFECT OF GENOTYPE, BLOCKING TX
y <- DGEList(counts = counts.keep.F)
y <- calcNormFactors(y)
Tx <- factor(sampleinfo.F$Tx)
Strain <- factor(sampleinfo.F$Strain)
design <- model.matrix(~Tx+Strain)

v <- voom(y,design,plot=TRUE)

fit <- lmFit(v, design)
fit_geno <- eBayes(fit)

Geno <- topTable(fit_geno,coef = "Strain1", n = "inf", sort = "p", adjust.method = "BH")

Geno$abs.logFC <- abs(Geno$logFC)
write.table(Geno, file="Geno-Tx_Female J ref.txt", row.names = TRUE, sep="\t")
#At this point I have had to manually add the header "Gene" to the outputed text file
#file in col 1

#Insert your directory here
Geno.F <- read.delim("C:/Users/jbeierle/Dropbox/lab/BALBc/BALBc QTL Mapping Fixed Genotyping/For GITHUB/Liver RNAseq/Geno-Tx_Female J ref.txt")


#-------DEG analysis MALE mice------
y <- DGEList(counts.keep.M)
names(y)   # See what slots are stored in y
y$samples   # Library size information is stored in the samples slot


#----EFFECT OF GENOTYPE, BLOCKING TX
y <- DGEList(counts = counts.keep.M)
y <- calcNormFactors(y)
Tx <- factor(sampleinfo.M$Tx)
Strain <- factor(sampleinfo.M$Strain)
design <- model.matrix(~Tx+Strain)

v <- voom(y,design,plot=TRUE)

fit <- lmFit(v, design)
fit_geno <- eBayes(fit)

Geno <- topTable(fit_geno,coef = "Strain1", n = "inf", sort = "p", adjust.method = "BH")

Geno$abs.logFC <- abs(Geno$logFC)
write.table(Geno, file="Geno-Tx_Male J ref.txt", row.names = TRUE, sep="\t")
#At this point I have had to manually add the header "Gene" to the outputed text file
#file in col 1

#Insert your directory here
Geno.M <- read.delim("C:/Users/jbeierle/Dropbox/lab/BALBc/BALBc QTL Mapping Fixed Genotyping/For GITHUB/Liver RNAseq/Geno-Tx_Male J ref.txt")


#----Extracting PK genes from GO terms----
#Extract Genes from GO terms
library(org.Mm.eg.db)
library(GO.db)
go_id = GOID( GOTERM[ Term(GOTERM) == "xenobiotic metabolic process"])
get(go_id, org.Mm.egGO2ALLEGS)
allegs = get(go_id, org.Mm.egGO2ALLEGS)
genes = unlist(mget(allegs,org.Mm.egSYMBOL))
genes.XMP <- genes


go_id = GOID( GOTERM[ Term(GOTERM) == "xenobiotic transport"])
get(go_id, org.Mm.egGO2ALLEGS)
allegs = get(go_id, org.Mm.egGO2ALLEGS)
genes = unlist(mget(allegs,org.Mm.egSYMBOL))
genes.Xt <- genes

genes.complete <- c(genes.XMP, genes.Xt)
genes.complete <- unique(genes.complete)
write_xlsx(as.data.frame(genes.complete), path="GOgenesOut.xlsx")

GOgenesOut <- read_excel("GOgenesOut.xlsx")
GO.genes <- GOgenesOut

#Removing genes added to PK gene list below
GO.genes <- GO.genes[str_detect(GO.genes$genes.complete, "Cyp", negate = TRUE), ]
GO.genes <- GO.genes[str_detect(GO.genes$genes.complete, "Abc", negate = TRUE), ]
GO.genes <- GO.genes[str_detect(GO.genes$genes.complete, "Ugt", negate = TRUE), ]
GO.genes <- GO.genes[str_detect(GO.genes$genes.complete, "Slc", negate = TRUE), ]


#----Sex Collapsed PK Genes----
Geno.All.Cyp <- Geno.All[str_detect(Geno.All$Gene, "Cyp"), ]
Geno.All.Abc <- Geno.All[str_detect(Geno.All$Gene, "Abc"), ]
Geno.All.Ugt <- Geno.All[str_detect(Geno.All$Gene, "Ugt"), ]
Geno.All.Slc <- Geno.All[str_detect(Geno.All$Gene, "Slc"), ]
Geno.All.GO <- c()


i = 1
for (i in 1:nrow(GO.genes)){
  temp <- GO.genes$genes.complete[i]
  j = 1
  
  for (j in 1:nrow(Geno.All)){
    if(isTRUE(temp == Geno.All$genes.complete[j])){
      Geno.All.GO <- rbind(Geno.All.GO, Geno.All[j,])
      break
    } else{NULL}}
  
}
Geno.All.PK <- rbind(Geno.All.Abc, Geno.All.Cyp, Geno.All.Slc, Geno.All.Ugt, Geno.All.GO)

#----Female Only PK Genes----
Geno.F.Cyp <- Geno.F[str_detect(Geno.F$Gene, "Cyp"), ]
Geno.F.Abc <- Geno.F[str_detect(Geno.F$Gene, "Abc"), ]
Geno.F.Ugt <- Geno.F[str_detect(Geno.F$Gene, "Ugt"), ]
Geno.F.Slc <- Geno.F[str_detect(Geno.F$Gene, "Slc"), ]
Geno.F.GO <- c()
i = 1
for (i in 1:nrow(GO.genes)){
  temp <- GO.genes$genes.complete[i]
  j = 1
  
  for (j in 1:nrow(Geno.F)){
    if(isTRUE(temp == Geno.F$genes.complete[j])){
      Geno.F.GO <- rbind(Geno.F.GO, Geno.F[j,])
      break
    } else{NULL}}
  
}
Geno.F.PK <- rbind(Geno.F.Abc, Geno.F.Cyp, Geno.F.Slc, Geno.F.Ugt, Geno.F.GO)


#----Male Only PK Genes----
Geno.M.Cyp <- Geno.M[str_detect(Geno.M$Gene, "Cyp"), ]
Geno.M.Abc <- Geno.M[str_detect(Geno.M$Gene, "Abc"), ]
Geno.M.Ugt <- Geno.M[str_detect(Geno.M$Gene, "Ugt"), ]
Geno.M.Slc <- Geno.M[str_detect(Geno.M$Gene, "Slc"), ]
Geno.M.GO <- c()
i = 1
for (i in 1:nrow(GO.genes)){
  temp <- GO.genes$genes.complete[i]
  j = 1
  
  for (j in 1:nrow(Geno.M)){
    if(isTRUE(temp == Geno.M$genes.complete[j])){
      Geno.M.GO <- rbind(Geno.M.GO, Geno.M[j,])
      break
    } else{NULL}}
  
}
Geno.M.PK <- rbind(Geno.M.Abc, Geno.M.Cyp, Geno.M.Slc, Geno.M.Ugt, Geno.M.GO)

#----JPET OMOR Figure 5----
#Figure 5A
#General DEG Plots
mutated_Geno.All <- mutate(Geno.All, sig=ifelse(abs(Geno.All$logFC)>.32, ifelse(Geno.All$P.Value<0.05, "both", "abs.logFC>0.32"), ifelse(Geno.All$P.Value<0.05,"p<0.05","not sig")))

volc = ggplot(mutated_Geno.All, aes(logFC, -log10(P.Value))) + 
  geom_point(aes(col=sig), size=2.5) + #add points colored by significance
  scale_color_manual(values=c("orange2","mediumseagreen","black", "royalblue3")) + 
  ggtitle("Differentially expressed genes within the Liver") +
  theme(legend.position="bottom", legend.title = element_blank(),
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14))
volc+geom_text_repel(data=subset(mutated_Geno.All, mutated_Geno.All$abs.logFC>0.5 | -log10(mutated_Geno.All$P.Value)>15), aes(label=Gene), size=4.5)


#Figure 5B
#Female PK genes
mutated_Geno.F.PK <- mutate(Geno.F.PK, sig=ifelse(abs(Geno.F.PK$logFC)>.32, ifelse(Geno.F.PK$P.Value<0.05, "both", "abs.logFC>0.32"), ifelse(Geno.F.PK$P.Value<0.05,"p<0.05","not sig")))

volc = ggplot(mutated_Geno.F.PK, aes(logFC, -log10(P.Value))) + 
  geom_point(aes(col=sig), size=2.5) + #add points colored by significance
  scale_color_manual(values=c("orange2","mediumseagreen","black", "royalblue3")) + 
  ggtitle("Female BALBc Liver PK Genes by Strain") +
  theme(legend.position="bottom", legend.title = element_blank(),
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14))
volc+geom_text_repel(data=subset(mutated_Geno.F.PK, mutated_Geno.F.PK$abs.logFC>0.5 | -log10(mutated_Geno.F.PK$P.Value)>15), aes(label=Gene), size=4.5)


#Figure 5C
#Male PK genes
mutated_Geno.M.PK <- mutate(Geno.M.PK, sig=ifelse(abs(Geno.M.PK$logFC)>.32, ifelse(Geno.M.PK$P.Value<0.05, "both", "abs.logFC>0.32"), ifelse(Geno.M.PK$P.Value<0.05,"p<0.05","not sig")))

volc = ggplot(mutated_Geno.M.PK, aes(logFC, -log10(P.Value))) + 
  geom_point(aes(col=sig), size=2.5) + #add points colored by significance
  scale_color_manual(values=c("orange2","mediumseagreen","black", "royalblue3")) + 
  ggtitle("Female BALBc Liver PK Genes by Strain") +
  theme(legend.position="bottom", legend.title = element_blank(),
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        legend.text = element_text(size=14))
volc+geom_text_repel(data=subset(mutated_Geno.M.PK, mutated_Geno.M.PK$abs.logFC>0.5 | -log10(mutated_Geno.M.PK$P.Value)>15), aes(label=Gene), size=4.5)