#----Load Packages----
library(splitstackshape)
library(tidyr)
library(MatrixEQTL)
library(reshape2)
library(XLConnect)
library(xlsx)
install.packages("XLConnect", type="source", INSTALL_opts = c("--no-multiarch"))
require(devtools)
install_version("XLConnect", version = "1.0.2", repos = "http://cran.us.r-project.org")
library(edgeR)
library(limma)
library(Glimma)
library(dplyr)
library(readxl)
library(calibrate)
library(HGNChelper)

#-----STR: Adjusting Striatial count file------
#Your directory here
count_str <- read.delim("C:/Users/jbeierle/Dropbox/lab/BALBc/BALBc QTL Mapping Fixed Genotyping/For GITHUB/eQTL Analysis/count_str.txt")
count_working<-count_str
#count_working<-count[,-c(8),]
colnames(count_working) <- c("Gene", "701", "702", "703", "704", "705", "706", "707", "708",
                             "712", "714", "715", "716", "718", "719", "721", "722", "723",
                             "725", "726", "727", "729", "733", "737", "738", "739", "744",
                             "745", "746", "747", "748", "801", "803", "805", "806", "807",
                             "808", "809", "810", "811", "812", "813", "815", "816", "817",
                             "818", "819", "821", "822", "825", "826", "827", "829", "830",
                             "831", "832", "833", "834", "835", "836", "838", "841", "842",
                             "844", "845")


#Remove unannotated genes
count_working<-count_working[!grepl("Gm", count_working$Gene),]
count_working<-count_working[!grepl("Rik", count_working$Gene),]

rownames(count_working) <- count_working[,1] #Make genes rownames
count_working <- count_working[,-(1)]

countdata.str <- count_working

#write.table(count_edited, file="count_edited_noGMs_Hcamp.txt", sep = "\t", quote=FALSE, row.name=FALSE, col.names = TRUE)




#Remove lowly expressed genes
myCPM <- cpm (countdata.str)
head(myCPM)
thresh <- myCPM > 0.30
head(thresh)
table(rowSums(thresh))
#keep genes that have at least 3 TRUES in each row of thresh - meaning genes that have >10 cpm in at least 3 samples
keep <- rowSums(thresh) >= 3
summary(keep)
# Subset the rows of countdata.str to keep the more highly expressed genes
counts.keep <- countdata.str[keep,]
dim(counts.keep)

write.table(counts.keep, file="counts_FinalQC_Str.txt", sep = "\t", quote=FALSE,
            row.name=TRUE, col.names = TRUE)

#-----STR: Setting Directories and parameters for analysis------
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelANOVA;

# Genotype, Expression, Covariate file name
SNP_file_name = paste( "geno_edited.txt", sep="");
snps_location_file_name = paste("marker_loc.txt", sep="");

expression_file_name = paste("counts_FinalQC_Str.txt", sep="");
gene_location_file_name = paste("gene_loc.txt", sep="");

covariates_file_name = paste("CovMartix.txt", sep="");  # Set to character() for no covariates

# Output file name
output_file_name = tempfile();
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-3;
pvOutputThreshold_cis = .001;
pvOutputThreshold_tra = .001;  ###if set to 0, turns off trans, if set to 1 turns off cutoff

# Error covariance matrix, Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


# Distance for local gene-SNP pairs
cisDist = 70e8; ###1.95e8 is the length of ch1


#-----STR: Loading Data------

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "-"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);


## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

#-----STR: Run the analysis------

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)
AllQTLs <- me$all$eqtls
#Export all eQTLs p<0.001
write.xlsx(me$all$eqtls, file = "All eqtls e-3 Str.xlsx",
sheetName = "all", append = FALSE)
write.xlsx(me$all$eqtls, file = "All eqtls e-3 Hcamp 10-7-2021.xlsx",
           sheetName = "all", append = FALSE)


#-----HIPP: Adjusting Striatial count file------
#Your directory here
count_hipp <- read.delim("C:/Users/jbeierle/Dropbox/lab/BALBc/BALBc QTL Mapping Fixed Genotyping/For GITHUB/eQTL Analysis/count_hipp.txt")
count_working<-count_hipp
#count_working<-count[,-c(8),]
colnames(count_working) <- c("Gene", "701", "702", "703", "704", "705", "706", "707", "708",
                             "712", "714", "715", "716", "718", "719", "721", "722", "723",
                             "725", "726", "727", "729", "733", "737", "738", "739", "744",
                             "745", "746", "747", "748", "801", "803", "805", "806", "807",
                             "808", "809", "810", "811", "812", "813", "815", "816", "817",
                             "818", "819", "821", "822", "825", "826", "827", "829", "830",
                             "831", "832", "833", "834", "835", "836", "838", "841", "842",
                             "844", "845")
#Remove unannotated genes
count_working<-count_working[!grepl("Gm", count_working$Gene),]
count_working<-count_working[!grepl("Rik", count_working$Gene),]
count_working<-count_working[!grepl("1-Mar", count_working$Gene),]
count_working<-count_working[!grepl("2-Mar", count_working$Gene),]

rownames(count_working) <- count_working[,1] #Make genes rownames
count_working <- count_working[,-(1)]

countdata.hipp <- count_working


#Remove lowly expressed genes
myCPM <- cpm (countdata.hipp)
head(myCPM)
thresh <- myCPM > 0.30
head(thresh)
table(rowSums(thresh))
#keep genes that have at least 3 TRUES in each row of thresh - meaning genes that have >10 cpm in at least 3 samples
keep <- rowSums(thresh) >= 3
summary(keep)
# Subset the rows of countdata.hipp to keep the more highly expressed genes
counts.keep <- countdata.hipp[keep,]
dim(counts.keep)

write.table(counts.keep, file="counts_FinalQC_Hipp.txt", sep = "\t", quote=FALSE,
            row.name=TRUE, col.names = TRUE)

#-----HIPP: Setting Directories and parameters for analysis------
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelANOVA;

# Genotype, Expression, Covariate file name
SNP_file_name = paste( "geno_edited.txt", sep="");
snps_location_file_name = paste("marker_loc.txt", sep="");

expression_file_name = paste("counts_FinalQC_Hipp.txt", sep="");
gene_location_file_name = paste("gene_loc.txt", sep="");

covariates_file_name = paste("CovMartix.txt", sep="");  # Set to character() for no covariates

# Output file name
output_file_name = tempfile();
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-3;
pvOutputThreshold_cis = .001;
pvOutputThreshold_tra = .001;  ###if set to 0, turns off trans, if set to 1 turns off cutoff

# Error covariance matrix, Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


# Distance for local gene-SNP pairs
cisDist = 70e8; ###1.95e8 is the length of ch1


#-----HIPP: Loading Data------

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "-"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);


## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

#-----STR: Run the analysis------

me = Matrix_eQTL_engine(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)
AllQTLs <- me$all$eqtls
#Export all eQTLs p<0.001

write.xlsx(me$all$eqtls, file = "All eqtls e-3 Hcamp.xlsx",
           sheetName = "all", append = FALSE)